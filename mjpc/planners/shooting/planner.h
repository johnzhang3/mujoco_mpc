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

#ifndef MJPC_PLANNERS_SHOOTING_PLANNER_H_
#define MJPC_PLANNERS_SHOOTING_PLANNER_H_

#include <vector>
#include <shared_mutex>

#include <mujoco/mujoco.h>
#include "mjpc/planners/planner.h"
#include "mjpc/planners/shooting/policy.h"
#include "mjpc/planners/shooting/gauss_newton.h"
#include "mjpc/planners/shooting/settings.h"
#include "mjpc/planners/model_derivatives.h"
#include "mjpc/planners/cost_derivatives.h"
#include "mjpc/states/state.h"
#include "mjpc/task.h"
#include "mjpc/threadpool.h"
#include "mjpc/trajectory.h"

namespace mjpc {

// Shooting-based planner using spline control parameterization
class ShootingPlanner : public Planner {
 public:
  ShootingPlanner() = default;
  ~ShootingPlanner() override = default;

  // Initialize data and settings
  void Initialize(mjModel* model, const Task& task) override;

  // Allocate memory
  void Allocate() override;

  // Reset memory to zeros
  void Reset(int horizon, const double* initial_repeated_action = nullptr) override;

  // Set state
  void SetState(const State& state) override;

  // Optimize policy using shooting method
  void OptimizePolicy(int horizon, ThreadPool& pool) override;

  // Compute trajectory using current policy
  void NominalTrajectory(int horizon, ThreadPool& pool) override;

  // Set action from policy
  void ActionFromPolicy(double* action, const double* state, double time, 
                       bool use_previous = false) override;

  // Return trajectory with best total return
  const Trajectory* BestTrajectory() override;

  // Visualize planner-specific traces
  void Traces(mjvScene* scn) override;

  // Planner-specific GUI elements  
  void GUI(mjUI& ui) override;

  // Planner-specific plots
  void Plots(mjvFigure* fig_planner, mjvFigure* fig_timer,
             int planner_shift, int timer_shift, int planning, int* shift) override;

  // Return number of parameters optimized by planner
  int NumParameters() override {
    return policy.spline_params.size();
  }

  // single gauss-newton iteration
  void Iteration(int horizon, ThreadPool& pool);

  // linesearch over action improvement
  void ActionRollouts(int horizon, ThreadPool& pool);

  // return index of trajectory with best rollout
  int BestRollout();

  void UpdateNumTrajectoriesFromGUI();

  // Members
  mjModel* model;
  const Task* task;

  // State
  std::vector<double> state;
  double time;
  std::vector<double> mocap;
  std::vector<double> userdata;

  // Policy
  ShootingPolicy policy;
  ShootingPolicy previous_policy;
  ShootingPolicy candidate_policy[kMaxTrajectory];
  
  // Dimensions
  int dim_state;
  int dim_state_derivative;
  int dim_action;
  int dim_sensor;
  int dim_max;

  // Trajectory
  Trajectory trajectory[kMaxTrajectory];

  // Derivatives
  ModelDerivatives model_derivative;
  CostDerivatives cost_derivative;

  // Optimization
  ShootingGaussNewton gauss_newton;

  // Settings
  ShootingSettings settings;

  // step sizes
  double linesearch_steps[kMaxTrajectory];

  // best trajectory id
  int winner;

  // values
  double action_step;
  double feedback_scaling;
  double improvement;
  double expected;
  double surprise;

  // Compute times
  double nominal_compute_time;
  double model_derivative_compute_time;
  double cost_derivative_compute_time;
  double optimization_compute_time;
  double rollouts_compute_time;
  double policy_update_compute_time;

  // Mutex for thread safety
  mutable std::shared_mutex mtx_;

 private:
  int num_trajectory_ = 1;
  int num_rollouts_gui_ = 1;

};

}  // namespace mjpc

#endif  // MJPC_PLANNERS_SHOOTING_PLANNER_H_