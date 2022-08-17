#ifndef INEKF_SLIP_ESTIMATOR_HPP_
#define INEKF_SLIP_ESTIMATOR_HPP_

#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/StdVector"

#include "legged_state_estimator/macros.hpp"
#include "legged_state_estimator/robot_model.hpp"
#include "legged_state_estimator/contact_estimator.hpp"
#include "legged_state_estimator/low_pass_filter.hpp"


namespace legged_state_estimator {

struct SlipEstimatorSettings {
  std::vector<double> beta0;
  std::vector<double> beta1;
  double slip_velocity_cov_alpha;
  double slip_prob_threshold;
  double lpf_contact_surface_normal_cutoff;
  double lpf_friction_coefficient_cutoff;
};


class SlipEstimator {
public:
  using Vector1d = Eigen::Matrix<double, 1, 1>;

  SlipEstimator(const RobotModel& robot_model, 
                const SlipEstimatorSettings& settings, const double dt);

  SlipEstimator();

  INEKF_USE_DEFAULT_DESTTUCTOR(SlipEstimator);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(SlipEstimator);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(SlipEstimator);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(SlipEstimator);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(SlipEstimator);

  void reset();

  // M(q) a + h (q, v) = S^T u + J^T f

  void update(const RobotModel& robot_model, 
              const ContactEstimator& contact_estimator);

  void setParameters(const SlipEstimatorSettings& settings);

  const std::vector<std::pair<int, bool>>& getSlipState() const;

  const std::vector<double>& getSlipProbability() const;

  const std::vector<double>& getSlipVelocityNorm() const;

  const std::vector<double>& getSlipVelocityCovariance() const;

  const std::vector<double>& getFrictionCoefficientEstimate() const;

  const std::vector<Eigen::Vector3d>& getContactSurfaceNormalEstimate() const;

  const std::vector<Eigen::Matrix3d>& getContactSurfaceEstimate() const;

  void resetContactSurfaceNormalEstimate(
      const std::vector<Eigen::Vector3d>& contact_surface_normal);

  void resetContactSurfaceNormalEstimate(
      const Eigen::Vector3d& contact_surface_normal=(Eigen::Vector3d() << 0., 0., 1.0).finished());

  void resetFrictionCoefficientEstimate(
      const std::vector<double>& friction_coefficient);

  void resetFrictionCoefficientEstimate(const double friction_coefficient=0.6);

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const SlipEstimator& d);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SlipEstimatorSettings settings_;
  std::vector<Eigen::Matrix3d> contact_surface_estimate_;
  std::vector<Eigen::Vector3d> contact_velocity_, contact_velocity_prev_, 
                               force_velocity_plane_normal_estimate_,
                               contact_surface_normal_estimate_,
                               contact_force_surface_local_;
  std::vector<double> contact_velocity_norm_, slip_probability_, slip_covariance_,
                      friction_coefficient_estimate_;
  std::vector<std::pair<int, bool>> slip_state_;
  std::vector<LowPassFilter<double, 3>> lpf_contact_surface_normal_;
  std::vector<LowPassFilter<double, 1>> lpf_friction_coefficient_;
  int num_contacts_;
};

} // namespace legged_state_estimator

#endif // INEKF_SLIP_ESTIMATOR_HPP_