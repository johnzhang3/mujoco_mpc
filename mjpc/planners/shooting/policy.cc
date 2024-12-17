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

#include "mjpc/planners/shooting/policy.h"
#include "mjpc/spline/spline.h"
#include <algorithm>
#include <mujoco/mujoco.h>
#include "mjpc/utilities.h"

namespace mjpc {

void ShootingPolicy::Allocate(const mjModel* model, const Task& task, int horizon) {
  // model
  this->model = model;

  // reference trajectory 
  trajectory.Initialize(model->nq + model->nv + model->na, model->nu,
                      task.num_residual, task.num_trace, kMaxTrajectoryHorizon);
  trajectory.Allocate(kMaxTrajectoryHorizon);

  // number of knots
  num_knots = GetNumberOrDefault(4, model, "shooting_num_knots");

  // spline parameters ((num_knots-1) * dim_action since first knot is fixed)
  spline_params.resize((num_knots - 1) * model->nu);

  // Initialize spline with model->nu dimensions
  plan = mjpc::spline::TimeSpline(model->nu);
  plan.Reserve(num_knots);

  // Set interpolation type from model or default
  plan.SetInterpolation(static_cast<mjpc::spline::SplineInterpolation>(
      GetNumberOrDefault(mjpc::spline::SplineInterpolation::kCubicSpline, 
                        model, "shooting_interpolation")));
}

void ShootingPolicy::Reset(int horizon, const double* initial_repeated_action) {
  // Reset trajectory including times
  trajectory.Reset(horizon, initial_repeated_action);
  
  // Set time points (equally spaced)
  for (int i = 0; i < horizon; i++) {
    trajectory.times[i] = i * model->opt.timestep;
  }
  
  // Reset spline parameters and plan
  if (initial_repeated_action) {
    // Set all knot points to initial action
    for (int i = 0; i < num_knots-1; i++) {
      mju_copy(spline_params.data() + i * model->nu, 
               initial_repeated_action, model->nu);
    }
  } else {
    std::fill(spline_params.begin(), spline_params.end(), 0.0);
  }

  // Reset spline plan
  plan.Clear();
  
  // Add initial control point
  plan.AddNode(0.0);
}

void ShootingPolicy::UpdateSplinePlan() const {
  // Clear existing nodes
  plan.Clear();
  
  double time_horizon = trajectory.times[trajectory.horizon-1];
  double dt = time_horizon / (num_knots - 1);

  // Add knot points to plan
  for (int i = 0; i < num_knots-1; i++) {
    auto node = plan.AddNode(i * dt);
    std::copy(spline_params.data() + i * model->nu, 
              spline_params.data() + (i + 1) * model->nu,
              node.values().begin());
  }
}

void ShootingPolicy::Action(double* action, const double* state, double time) const {
  // Get time horizon
  double time_horizon = trajectory.times[trajectory.horizon-1];
  if (time_horizon <= 0) {
    mju_zero(action, model->nu);
    return;
  }

  // For first timestep, use first knot point directly
  if (time <= model->opt.timestep) {
    mju_copy(action, spline_params.data(), model->nu);
    Clamp(action, model->actuator_ctrlrange, model->nu);
    return;
  }

  // Update spline plan with current parameters
  UpdateSplinePlan();

  // Sample spline at current time
  plan.Sample(time, absl::MakeSpan(action, model->nu));

  // Clamp controls
  Clamp(action, model->actuator_ctrlrange, model->nu);
}

void ShootingPolicy::CopyFrom(const ShootingPolicy& policy, int horizon) {
  // Copy trajectory
  trajectory = policy.trajectory;

  // Copy spline parameters
  mju_copy(spline_params.data(), policy.spline_params.data(), 
           (num_knots-1) * model->nu);

  // Copy spline plan configuration
  plan = policy.plan;
}
}  // namespace mjpc