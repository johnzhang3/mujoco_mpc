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

#ifndef MJPC_PLANNERS_SHOOTING_GAUSS_NEWTON_H_
#define MJPC_PLANNERS_SHOOTING_GAUSS_NEWTON_H_

#include <vector>
#include <memory>
#include "mjpc/planners/shooting/settings.h"
#include "mjpc/planners/cost_derivatives.h"
#include "mjpc/planners/model_derivatives.h"
#include "mjpc/planners/shooting/policy.h"
#include "mjpc/task.h"

namespace mjpc {

// Data and methods for Gauss-Newton optimization step
class ShootingGaussNewton {
 public:
  ShootingGaussNewton() = default;
  ~ShootingGaussNewton() = default;

  // Allocate memory
  void Allocate(int dim_state_derivative, int dim_action, int num_knots);

  // Reset memory to zeros
  void Reset(int dim_state_derivative, int dim_action, int num_knots);

  // Compute gradient and Gauss-Newton step
  int OptimizationStep(const ModelDerivatives* md, const CostDerivatives* cd,
                      int dim_state_derivative, int dim_action,
                      double* params, double lambda, int num_knots,
                      const double* actions, const double* action_limits,
                      const ShootingSettings& settings);

  // Scale regularization
  void ScaleRegularization(double factor, double reg_min, double reg_max);

  // Update regularization based on improvement ratio
  void UpdateRegularization(double reg_min, double reg_max, double z, double s);

  // Compute rollout cost for given parameters
  double RolloutCost(const std::vector<double>& params,
                    const double* state,
                    double time,
                    const double* mocap,
                    const double* userdata,
                    int horizon);

  // Required references for optimization
  ShootingPolicy* policy;        // Policy being optimized
  const mjModel* model;          // MuJoCo model
  const Task* task;              // Task specification
  std::vector<std::unique_ptr<mjData>> data;  // MuJoCo data for rollouts
  std::vector<double> mocap;     // Mocap data for rollouts
  std::vector<double> userdata;  // User data for rollouts

  // Members for optimization
  std::vector<double> gradient;      // Gradient of cost wrt parameters
  std::vector<double> hessian;       // Gauss-Newton approximate Hessian
  std::vector<double> step;          // Optimization step
  std::vector<double> params_scratch;  // Parameter perturbation scratch space
  std::vector<double> grad_scratch;    // Gradient computation scratch space
  double dV[2];                      // Expected vs actual cost reduction

  // Regularization parameters  
  double regularization;          // Current regularization value
  double regularization_rate;     // Rate of regularization change  
  double regularization_factor;   // Factor for scaling regularization
};

}  // namespace mjpc

#endif  // MJPC_PLANNERS_SHOOTING_GAUSS_NEWTON_H_