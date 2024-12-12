#ifndef MJPC_PLANNERS_SHOOTING_SETTINGS_H_
#define MJPC_PLANNERS_SHOOTING_SETTINGS_H_

namespace mjpc {

// Settings for Shooting-based Gauss-Newton optimization
struct ShootingSettings {
  double fd_tolerance = 1.0e-6;   // finite difference tolerance for gradient
  double fd_mode = 0;  // type of finite difference; 0: one-sided, 1: centered
  double min_step_size = 1.0e-3;  // minimum step size (line search)
  double lambda_init = 1.0e-3;    // initial regularization factor
  double lambda_max = 1.0e6;      // maximum regularization
  double lambda_factor = 2.0;     // factor by which regularization is scaled up/down
  int max_iterations = 50;        // maximum Gauss-Newton iterations
  int verbose = 1;                // print debug info

  // If control constraints are needed, add similar parameters as in iLQG.
  int action_limits = 1; // whether to enforce actuator_ctrlrange
};

}  // namespace mjpc

#endif  // MJPC_PLANNERS_SHOOTING_SETTINGS_H_
