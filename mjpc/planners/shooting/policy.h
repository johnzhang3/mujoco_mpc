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

#ifndef MJPC_PLANNERS_SHOOTING_POLICY_H_
#define MJPC_PLANNERS_SHOOTING_POLICY_H_

#include <vector>
#include <mujoco/mujoco.h>
#include "mjpc/planners/policy.h"
#include "mjpc/task.h"
#include "mjpc/trajectory.h"
#include "mjpc/spline/spline.h" 

namespace mjpc {

// Shooting policy using spline parameterization
class ShootingPolicy : public Policy {
 public:
  ShootingPolicy() = default;
  ~ShootingPolicy() override = default;

  // Allocate memory
  void Allocate(const mjModel* model, const Task& task, int horizon) override;

  // Reset memory to zeros
  void Reset(int horizon, const double* initial_repeated_action = nullptr) override;

  // Set action from policy
  void Action(double* action, const double* state, double time) const override;

  // Copy policy
  void CopyFrom(const ShootingPolicy& policy, int horizon);

  // Update spline plan with current parameters
  void UpdateSplinePlan() const;

  // Members
  const mjModel* model;
  Trajectory trajectory;                // Reference trajectory
  std::vector<double> spline_params;   // Spline control points
  std::vector<double> control_scratch; // Scratch space for interpolation
  int num_knots;                      // Number of spline control points
  mutable spline::TimeSpline plan;     // Added spline plan
};

}  // namespace mjpc

#endif  // MJPC_PLANNERS_SHOOTING_POLICY_H_