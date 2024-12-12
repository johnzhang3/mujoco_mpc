// Copyright 2024 DeepMind Technologies Limited
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "mjpc/planners/shooting/planner.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <iostream>

#include <mujoco/mujoco.h>
#include "mjpc/array_safety.h"
#include "mjpc/utilities.h"

namespace mjpc {
namespace mju = ::mujoco::util_mjpc;

// Initialize data and settings
void ShootingPlanner::Initialize(mjModel* model, const Task& task) {
  // delete mjData instances since model might have changed
  data_.clear();
  
  // allocate one mjData for nominal
  ResizeMjData(model, 1);
  
  // model
  this->model = model;
  
  // task
  this->task = &task;
  
  // dimensions
  dim_state = model->nq + model->nv + model->na;
  dim_state_derivative = 2 * model->nv + model->na;
  dim_action = model->nu;
  dim_sensor = model->nsensordata;
  dim_max = mju_max(mju_max(mju_max(dim_state, dim_state_derivative),
                           dim_action), model->nuser_sensor);

  // Read time step from model or use default
  double time_step = GetNumberOrDefault(0.1, model, "shooting_time_step");
  model->opt.timestep = time_step;
}

// Allocate memory
void ShootingPlanner::Allocate() {
  // state
  state.resize(model->nq + model->nv + model->na);
  mocap.resize(7 * model->nmocap);
  userdata.resize(model->nuserdata);
  
  // trajectories 
  for (int i = 0; i < kMaxTrajectory; i++) {
    trajectory[i].Initialize(dim_state, dim_action, task->num_residual,
                          task->num_trace, kMaxTrajectoryHorizon);
    trajectory[i].Allocate(kMaxTrajectoryHorizon);
  }
  
  // model derivatives
  model_derivative.Allocate(dim_state_derivative, dim_action,
                           dim_sensor, kMaxTrajectoryHorizon);
  
  // cost derivatives
  cost_derivative.Allocate(dim_state_derivative, dim_action,
                          task->num_residual, kMaxTrajectoryHorizon, dim_max);
  
  // gauss-newton step
  gauss_newton.Allocate(dim_state_derivative, dim_action,
                        GetNumberOrDefault(4, model, "shooting_num_knots"));
  
  // policy
  policy.Allocate(model, *task, kMaxTrajectoryHorizon);
  previous_policy.Allocate(model, *task, kMaxTrajectoryHorizon);
  for (int i = 0; i < kMaxTrajectory; i++) {
    candidate_policy[i].Allocate(model, *task, kMaxTrajectoryHorizon);
  }
}

// Reset memory to zeros
void ShootingPlanner::Reset(int horizon, const double* initial_repeated_action) {
  // state
  std::fill(state.begin(), state.end(), 0.0);
  std::fill(mocap.begin(), mocap.end(), 0.0);
  std::fill(userdata.begin(), userdata.end(), 0.0);
  time = 0.0;
  
  // trajectories
  for (int i = 0; i < kMaxTrajectory; i++) {
    trajectory[i].Reset(horizon, initial_repeated_action);
  }
  
  // model derivatives
  model_derivative.Reset(dim_state_derivative, dim_action,
                        dim_sensor, horizon);
  
  // cost derivatives
  cost_derivative.Reset(dim_state_derivative, dim_action,
                       task->num_residual, horizon);
  
  // gauss-newton step
  gauss_newton.Reset(dim_state_derivative, dim_action, policy.num_knots);
  
  // policies
  policy.Reset(horizon, initial_repeated_action);
  previous_policy.Reset(horizon, initial_repeated_action);
  for (int i = 0; i < kMaxTrajectory; i++) {
    candidate_policy[i].Reset(horizon, initial_repeated_action);
  }
}

// Set state
void ShootingPlanner::SetState(const State& state) {
  state.CopyTo(this->state.data(), this->mocap.data(),
               this->userdata.data(), &this->time);
}

// Optimize policy using shooting method
void ShootingPlanner::OptimizePolicy(int horizon, ThreadPool& pool) {
  // get nominal trajectory
  this->NominalTrajectory(horizon, pool);
  
  // iteration
  this->Iteration(horizon, pool);
}

// Compute trajectory using nominal policy
void ShootingPlanner::NominalTrajectory(int horizon, ThreadPool& pool) {
  if (num_trajectory_ == 0) {
    return;
  }
  // set policy trajectory horizon
  policy.trajectory.horizon = horizon;
  
  // resize data for rollouts
  ResizeMjData(model, pool.NumThreads());
  
  // start timer
  auto nominal_start = std::chrono::steady_clock::now();
  
  // rollout nominal trajectory
  {
    const std::shared_lock<std::shared_mutex> lock(mtx_);
    auto nominal_policy = [this](double* action, const double* state,
                                double time) {
      policy.Action(action, state, time);
    };
    
    policy.trajectory.Rollout(nominal_policy, task, model,
                             data_[0].get(), state.data(), time,
                             mocap.data(), userdata.data(), horizon);
  }
  
  // end timer
  nominal_compute_time = GetDuration(nominal_start);
}

void ShootingPlanner::UpdateNumTrajectoriesFromGUI() {
  num_trajectory_ = mju_min(num_rollouts_gui_, kMaxTrajectory);
}


void ShootingPlanner::Iteration(int horizon, ThreadPool& pool) {
  // previous cost
  double previous_cost = policy.trajectory.total_return;

  // resize data for rollouts
  ResizeMjData(model, pool.NumThreads());

  // ----- model derivatives ----- //
  auto model_derivative_start = std::chrono::steady_clock::now();

  model_derivative.Compute(model, data_, policy.trajectory.states.data(),
                         policy.trajectory.actions.data(),
                         policy.trajectory.times.data(), dim_state,
                         dim_state_derivative, dim_action, dim_sensor,
                         horizon, settings.fd_tolerance, settings.fd_mode,
                         pool, 0);

  model_derivative_compute_time = GetDuration(model_derivative_start);

  // ----- cost derivatives ----- //
  auto cost_derivative_start = std::chrono::steady_clock::now();

  cost_derivative.Compute(policy.trajectory.residual.data(),
                        model_derivative.C.data(),
                        model_derivative.D.data(),
                        dim_state_derivative, dim_action, dim_max,
                        dim_sensor, task->num_residual,
                        task->dim_norm_residual.data(),
                        task->num_term, task->weight.data(),
                        task->norm.data(), task->norm_parameter.data(),
                        task->num_norm_parameter.data(),
                        task->risk, horizon, pool);

  cost_derivative_compute_time = GetDuration(cost_derivative_start);

  // ----- optimization step ----- //
  auto optimization_start = std::chrono::steady_clock::now();

  int status = gauss_newton.OptimizationStep(
      &model_derivative, &cost_derivative,
      dim_state_derivative, dim_action,
      policy.spline_params.data(),
      gauss_newton.regularization,
      policy.num_knots,
      policy.trajectory.actions.data(),
      model->actuator_ctrlrange,
      settings);

  optimization_compute_time = GetDuration(optimization_start);

  if (!status) {
    if (settings.verbose) {
      std::cout << "Gauss-Newton step failed\n";
      std::cout << "  regularization: " << gauss_newton.regularization << "\n";
    }
    return;
  }

  // ----- rollout policy ----- //
  auto rollouts_start = std::chrono::steady_clock::now();

  // copy policy
  for (int j = 1; j < num_trajectory_; j++) {
    candidate_policy[j].CopyFrom(candidate_policy[0], horizon);
  }

  // ----- line search ----- //
  this -> ActionRollouts(horizon, pool);

  // get best rollout
  int best_rollout = this->BestRollout();
  std::cout << "Best rollout: " << best_rollout << std::endl;
  if (best_rollout == -1) {
    return;
  } else {
    winner = best_rollout;
  }

  // update nominal with winner
  candidate_policy[0].trajectory = trajectory[winner];
  // improvement
  action_step = linesearch_steps[winner];
  expected = -1.0 * action_step *
                  (gauss_newton.dV[0] + action_step * gauss_newton.dV[1]) +
              1.0e-16;
  improvement = previous_cost - trajectory[winner].total_return;
  surprise = mju_min(mju_max(0, improvement / expected), 2);

//   double actual_reduction = previous_cost - policy.trajectory.total_return;
//   double predicted_reduction = -gauss_newton.dV[0];
//   double reduction_ratio = improvement / (expected + 1e-6);

  rollouts_compute_time = GetDuration(rollouts_start);

    // update regularization
  gauss_newton.UpdateRegularization(settings.lambda_init,
                                  settings.lambda_max,
                                  surprise, action_step);

  if (settings.verbose) {
    std::cout << "\nShooting Iteration Info:\n";
    std::cout << "  best return: " << policy.trajectory.total_return << "\n";
    std::cout << "  previous return: " << previous_cost << "\n";
    std::cout << "  regularization: " << gauss_newton.regularization << "\n";
    std::cout << "  winner step size: " << linesearch_steps[winner] << "\n\n";
    std::cout << "Compute times (ms):\n";
    std::cout << "  model derivative: " << model_derivative_compute_time * 1e-3 << "\n";
    std::cout << "  cost derivative: " << cost_derivative_compute_time * 1e-3 << "\n";
    std::cout << "  rollouts: " << rollouts_compute_time * 1e-3 << "\n";
    std::cout << "  optimization: " << optimization_compute_time * 1e-3 << "\n\n";
  }
  // ----- policy update ----- //
  auto policy_update_start = std::chrono::steady_clock::now();

  {
    const std::shared_lock<std::shared_mutex> lock(mtx_);
    previous_policy = policy;
    policy.CopyFrom(candidate_policy[winner], horizon);
    
  }
  // stop timer
  policy_update_compute_time = GetDuration(policy_update_start);
}

// Adding ActionRollouts implementation to planner.cc:
void ShootingPlanner::ActionRollouts(int horizon, ThreadPool& pool) {
  // compute line search step sizes (log scaling)
  LogScale(linesearch_steps, 1.0, settings.min_step_size, num_trajectory_ - 1);
  linesearch_steps[num_trajectory_ - 1] = 0.0;  // One rollout with zero step

  // no one else should be writing
  {
    const std::shared_lock<std::shared_mutex> lock(mtx_);
    for (int i = 0; i < num_trajectory_; i++) {
      candidate_policy[i].CopyFrom(policy, horizon);
      
      // update spline parameters using Gauss-Newton step
      int num_params = candidate_policy[i].spline_params.size();
      for (int j = 0; j < num_params; j++) {
        candidate_policy[i].spline_params[j] += linesearch_steps[i] * gauss_newton.step[j];
      }

      // enforce control limits if enabled
      if (settings.action_limits) {
        for (int j = 0; j < candidate_policy[i].num_knots-1; j++) {
          int idx = j * model->nu;
          Clamp(candidate_policy[i].spline_params.data() + idx,
                model->actuator_ctrlrange,
                model->nu);
        }
      }
    }
  }

  // parallel rollouts
  int count_before = pool.GetCount();
  for (int i = 0; i < num_trajectory_; i++) {
    pool.Schedule([&data = data_, &trajectory = trajectory,
                  &candidate_policy = candidate_policy, &model = this->model,
                  &task = this->task, &state = this->state, &time = this->time,
                  &mocap = this->mocap, &userdata = this->userdata, 
                  horizon, i]() {  // Added horizon to capture list
      // policy rollout with candidate policy
      auto rollout_policy = [&candidate_policy = candidate_policy[i]](
                               double* action, const double* state, 
                               double time) {
        candidate_policy.Action(action, state, time);
      };

      // rollout trajectory using discrete time
      trajectory[i].RolloutDiscrete(rollout_policy, task, model,
                                  data[ThreadPool::WorkerId()].get(),
                                  state.data(), time, mocap.data(),
                                  userdata.data(), horizon);
    });
  }
  pool.WaitCount(count_before + num_trajectory_);
  pool.ResetCount();
}

int ShootingPlanner::BestRollout() {
  double best_return = 0;
  int best_rollout = -1;

  for (int j = num_trajectory_ - 1; j >= 0; j--) {
    if (trajectory[j].failure) continue;
    double rollout_return = trajectory[j].total_return;
    if (best_rollout == -1 || rollout_return < best_return) {
      best_return = rollout_return;
      best_rollout = j;
    }
  }

  if (settings.verbose) {
    for (int j = 0; j < num_trajectory_; j++) {
      std::cout << "Rollout " << j << " step: " << linesearch_steps[j] 
                << " return: " << trajectory[j].total_return 
                << (trajectory[j].failure ? " (failed)" : "") << std::endl;
    }
    if (best_rollout >= 0) {
      std::cout << "Best rollout: " << best_rollout 
                << " (return: " << best_return << ")" << std::endl;
    }
  }

  return best_rollout;
}

// Set action from policy
void ShootingPlanner::ActionFromPolicy(double* action, const double* state,
                                     double time, bool use_previous) {
  const std::shared_lock<std::shared_mutex> lock(mtx_);
  if (use_previous) {
    previous_policy.Action(action, state, time);
  } else {
    policy.Action(action, state, time);
  }
}

// Return trajectory with best total return
const Trajectory* ShootingPlanner::BestTrajectory() {
  const std::shared_lock<std::shared_mutex> lock(mtx_);
  return &policy.trajectory;
}


// Visualize planner-specific traces
void ShootingPlanner::Traces(mjvScene* scn) {
  // No special visualization needed
}

// GUI elements
void ShootingPlanner::GUI(mjUI& ui) {
  mjuiDef defShooting[] = {
    {mjITEM_SELECT, "Shooting Opt", 2, &settings.verbose, "Off\nOn"},
    {mjITEM_END}
  };
  mjui_add(&ui, defShooting);
}

// Planner plots
// Planner plots
void ShootingPlanner::Plots(mjvFigure* fig_planner, mjvFigure* fig_timer,
                           int planner_shift, int timer_shift, int planning,
                           int* shift) {
  // bounds
  double planner_bounds[2] = {-6, 6};

  // ----- planner ----- //

  // regularization
  mjpc::PlotUpdateData(fig_planner, planner_bounds,
                       fig_planner->linedata[0 + planner_shift][0] + 1,
                       mju_log10(mju_max(gauss_newton.regularization, 1.0e-6)),
                       100, 0 + planner_shift, 0, 1, -100);

  // best step size from line search
  double step_size = winner >= 0 ? linesearch_steps[winner] : 1.0e-6;
  mjpc::PlotUpdateData(fig_planner, planner_bounds,
                       fig_planner->linedata[1 + planner_shift][0] + 1,
                       mju_log10(mju_max(step_size, 1.0e-6)), 100,
                       1 + planner_shift, 0, 1, -100);

  // improvement ratio
  double actual_reduction = winner >= 0 ? 
      policy.trajectory.total_return - trajectory[winner].total_return : 0.0;
  double predicted_reduction = -gauss_newton.dV[0];
  double reduction_ratio = actual_reduction / (predicted_reduction + 1.0e-6);
  mjpc::PlotUpdateData(fig_planner, planner_bounds,
                       fig_planner->linedata[2 + planner_shift][0] + 1,
                       mju_log10(mju_max(reduction_ratio, 1.0e-6)), 100,
                       2 + planner_shift, 0, 1, -100);

  // legend
  mju::strcpy_arr(fig_planner->linename[0 + planner_shift], "Regularization");
  mju::strcpy_arr(fig_planner->linename[1 + planner_shift], "Step Size");
  mju::strcpy_arr(fig_planner->linename[2 + planner_shift], "Reduction Ratio");

  // ranges
  fig_planner->range[1][0] = planner_bounds[0];
  fig_planner->range[1][1] = planner_bounds[1];

  // ----- timer ----- //
  double timer_bounds[2] = {0, 1};

  // update plots
  PlotUpdateData(fig_timer, timer_bounds,
                 fig_timer->linedata[0 + timer_shift][0] + 1,
                 1.0e-3 * nominal_compute_time * planning, 100, 0 + timer_shift,
                 0, 1, -100);

  PlotUpdateData(fig_timer, timer_bounds,
                 fig_timer->linedata[1 + timer_shift][0] + 1,
                 1.0e-3 * model_derivative_compute_time * planning, 100,
                 1 + timer_shift, 0, 1, -100);

  PlotUpdateData(fig_timer, timer_bounds,
                 fig_timer->linedata[2 + timer_shift][0] + 1,
                 1.0e-3 * cost_derivative_compute_time * planning, 100,
                 2 + timer_shift, 0, 1, -100);

  PlotUpdateData(fig_timer, timer_bounds,
                 fig_timer->linedata[3 + timer_shift][0] + 1,
                 1.0e-3 * optimization_compute_time * planning, 100,
                 3 + timer_shift, 0, 1, -100);

  PlotUpdateData(fig_timer, timer_bounds,
                 fig_timer->linedata[4 + timer_shift][0] + 1,
                 1.0e-3 * rollouts_compute_time * planning, 100,
                 4 + timer_shift, 0, 1, -100);

  PlotUpdateData(fig_timer, timer_bounds,
                 fig_timer->linedata[5 + timer_shift][0] + 1,
                 1.0e-3 * policy_update_compute_time * planning, 100,
                 5 + timer_shift, 0, 1, -100);

  // legend
  mju::strcpy_arr(fig_timer->linename[0 + timer_shift], "Nominal");
  mju::strcpy_arr(fig_timer->linename[1 + timer_shift], "Model Deriv.");
  mju::strcpy_arr(fig_timer->linename[2 + timer_shift], "Cost Deriv.");
  mju::strcpy_arr(fig_timer->linename[3 + timer_shift], "Gauss-Newton");
  mju::strcpy_arr(fig_timer->linename[4 + timer_shift], "Rollouts");
  mju::strcpy_arr(fig_timer->linename[5 + timer_shift], "Policy Update");

  fig_timer->range[0][0] = -100;
  fig_timer->range[0][1] = 0;
  fig_timer->range[1][0] = 0.0;
  fig_timer->range[1][1] = timer_bounds[1];

  // planner shift
  shift[0] += 3;

  // timer shift
  shift[1] += 6;
}
}  // namespace mjpc