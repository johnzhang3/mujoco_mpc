#include "mjpc/planners/shooting/gauss_newton.h"
#include "mjpc/spline/spline.h"
#include <algorithm>
#include <mujoco/mujoco.h>
#include "mjpc/utilities.h"

namespace mjpc {

void ShootingGaussNewton::Allocate(int dim_state_derivative, int dim_action, 
                                  int num_knots) {
  // optimization variables - each knot point has dim_action parameters
  int num_params = (num_knots - 1) * dim_action;
  gradient.resize(num_params);
  hessian.resize(num_params * num_params);
  step.resize(num_params);
  
  // scratch space for finite differences
  params_scratch.resize(num_params);
  grad_scratch.resize(num_params);
  
  // reset values
  regularization = 1.0;
  regularization_rate = 1.0;
  regularization_factor = 2.0;
}

void ShootingGaussNewton::Reset(int dim_state_derivative, int dim_action, 
                               int num_knots) {
  std::fill(gradient.begin(), gradient.end(), 0.0);
  std::fill(hessian.begin(), hessian.end(), 0.0);
  std::fill(step.begin(), step.end(), 0.0);
  std::fill(params_scratch.begin(), params_scratch.end(), 0.0);
  std::fill(grad_scratch.begin(), grad_scratch.end(), 0.0);
  mju_zero(dV, 2);
  
  regularization = 1.0;
  regularization_rate = 1.0;
  regularization_factor = 2.0;
}

double ShootingGaussNewton::RolloutCost(
    const std::vector<double>& params,
    const double* state,
    double time,
    const double* mocap,
    const double* userdata,
    int horizon) {
    
  // Copy parameters to policy
  mju_copy(policy->spline_params.data(), params.data(), params.size());

  // Create policy callback
  auto rollout_policy = [this](double* action, const double* state, double time) {
    policy->Action(action, state, time);
  };

  // Rollout trajectory
  policy->trajectory.RolloutDiscrete(
      rollout_policy, task, model, data[0].get(),
      state, time, mocap, userdata, horizon);

  return policy->trajectory.total_return;
}
int ShootingGaussNewton::OptimizationStep(const ModelDerivatives* /*md*/,
    const CostDerivatives* /*cd*/, int dim_state_derivative, int dim_action,
    double* params, double lambda, int num_knots, const double* actions,
    const double* action_limits, const ShootingSettings& settings) {
    
  int num_params = (num_knots - 1) * dim_action;
  double eps = settings.fd_tolerance;

  // Get current state and other data from policy's trajectory for the rollout starting point
  const double* state = policy->trajectory.states.data();
  double time = policy->trajectory.times[0];
  int horizon = policy->trajectory.horizon;

  // Current parameters and nominal cost
  std::vector<double> params_current(params, params + num_params);
  
  // Create rollout policy callback
  auto rollout_policy = [this](double* action, const double* state, double time) {
    policy->Action(action, state, time);
  };

  // Compute nominal trajectory cost
  mju_copy(policy->spline_params.data(), params_current.data(), num_params);
  policy->trajectory.RolloutDiscrete(rollout_policy, task, model, 
                                   data[0].get(), state, time, 
                                   mocap.data(), userdata.data(), horizon);
  double nominal_cost = policy->trajectory.total_return;

  // Compute gradient with one-sided finite differences
  for (int i = 0; i < num_params; i++) {
    // Perturb parameter
    std::vector<double> params_perturbed = params_current;
    params_perturbed[i] += eps;
    
    // Copy perturbed parameters to policy
    mju_copy(policy->spline_params.data(), params_perturbed.data(), num_params);
    
    // Rollout perturbed trajectory
    policy->trajectory.RolloutDiscrete(rollout_policy, task, model,
                                     data[0].get(), state, time,
                                     mocap.data(), userdata.data(), horizon);
    double perturbed_cost = policy->trajectory.total_return;
    
    // Finite difference for gradient
    gradient[i] = (perturbed_cost - nominal_cost) / eps;
  }

  // Compute Gauss-Newton Hessian approximation
  for (int i = 0; i < num_params; i++) {
    std::vector<double> params_perturbed = params_current;
    params_perturbed[i] += eps;
    mju_copy(policy->spline_params.data(), params_perturbed.data(), num_params);
    
    // Get perturbed gradient column
    for (int j = 0; j < num_params; j++) {
      std::vector<double> params_double_perturbed = params_perturbed;
      params_double_perturbed[j] += eps;
      
      // Rollout double-perturbed trajectory
      mju_copy(policy->spline_params.data(), params_double_perturbed.data(), num_params);
      policy->trajectory.RolloutDiscrete(rollout_policy, task, model,
                                       data[0].get(), state, time,
                                       mocap.data(), userdata.data(), horizon);
      double cost_ij = policy->trajectory.total_return;
      
      // Rollout single-perturbed trajectory
      mju_copy(policy->spline_params.data(), params_perturbed.data(), num_params);
      policy->trajectory.RolloutDiscrete(rollout_policy, task, model,
                                       data[0].get(), state, time,
                                       mocap.data(), userdata.data(), horizon);
      double cost_i = policy->trajectory.total_return;
      
      // Second order finite difference
      grad_scratch[j] = (cost_ij - cost_i) / eps;
    }

    // Fill Hessian column
    for (int j = 0; j < num_params; j++) {
      hessian[j * num_params + i] = (grad_scratch[j] - gradient[j]) / eps;
    }
  }

  // Symmetrize Hessian
  for (int i = 0; i < num_params; i++) {
    for (int j = 0; j < i; j++) {
      double avg = 0.5 * (hessian[i * num_params + j] + hessian[j * num_params + i]);
      hessian[i * num_params + j] = hessian[j * num_params + i] = avg;
    }
  }

  // Add regularization to diagonal
  for (int i = 0; i < num_params; i++) {
    hessian[i * num_params + i] += lambda;
  }

  // Solve for step: H * step = -gradient
  mju_copy(step.data(), gradient.data(), num_params);
  if (mju_cholFactor(hessian.data(), num_params, 0.0) > 0) {
    mju_cholSolve(step.data(), hessian.data(), gradient.data(), num_params);
    mju_scl(step.data(), step.data(), -1.0, num_params);

    // Store predicted improvement 
    std::vector<double> temp(num_params);
    
    // Compute dV[0] = -gradient^T * step
    dV[0] = 0;
    for (int i = 0; i < num_params; i++) {
        dV[0] -= gradient[i] * step[i];
    }
    
    // Compute dV[1] = 0.5 * step^T * H * step
    // First compute H * step -> temp
    for (int i = 0; i < num_params; i++) {
        temp[i] = 0;
        for (int j = 0; j < num_params; j++) {
            temp[i] += hessian[i * num_params + j] * step[j];
        }
    }
    
    // Then compute step^T * temp
    dV[1] = 0;
    for (int i = 0; i < num_params; i++) {
        dV[1] += 0.5 * step[i] * temp[i];
    }

    // Restore nominal parameters
    mju_copy(policy->spline_params.data(), params_current.data(), num_params);
    return 1;
  }

  // Restore nominal parameters on failure
  mju_copy(policy->spline_params.data(), params_current.data(), num_params);
  return 0;
}

void ShootingGaussNewton::ScaleRegularization(double factor, double reg_min, 
                                             double reg_max) {
  if (factor > 1) {
    regularization_rate = mju_max(regularization_rate * factor, factor);
  } else {
    regularization_rate = mju_min(regularization_rate * factor, factor);
  }
  
  regularization = mju_min(mju_max(regularization * regularization_rate, 
                                  reg_min), reg_max);
}

void ShootingGaussNewton::UpdateRegularization(double reg_min, double reg_max,
                                             double z, double s) {
  // divergence or no improvement: increase regularization by factor^2  
  if (mju_isBad(z) || mju_isBad(s)) {
    this->ScaleRegularization(regularization_factor * regularization_factor,
                             reg_min, reg_max);
  } 
  // sufficient improvement: decrease regularization by factor
  else if (z > 0.5 || s > 0.3) {
    this->ScaleRegularization(1.0 / regularization_factor, reg_min, reg_max);
  }
  // insufficient improvement: increase regularization by factor
  else if (z < 0.1 || s < 0.06) {
    this->ScaleRegularization(regularization_factor, reg_min, reg_max);
  }
}

}  // namespace mjpc