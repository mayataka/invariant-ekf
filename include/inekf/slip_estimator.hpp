#ifndef INEKF_SLIP_ESTIMATOR_HPP_
#define INEKF_SLIP_ESTIMATOR_HPP_

#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/StdVector"

#include "inekf/macros.hpp"
#include "inekf/robot_model.hpp"
#include "inekf/contact_estimator.hpp"


namespace inekf {

struct SlipEstimatorSettings {
  std::vector<double> beta0;
  std::vector<double> beta1;
  double slip_velocity_cov_alpha;
  double slip_prob_threshold;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


class SlipEstimator {
public:
  SlipEstimator(const RobotModel& robot_model, 
                const SlipEstimatorSettings& settings);

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

  const std::vector<double>& getSlipVelocityCovariance() const;

  const std::vector<double>& getFrictionCoefficientEstimate() const;

  const std::vector<Eigen::Vector3d>& getContactSurfaceNormalEstimate() const;

  const std::vector<Eigen::Matrix3d>& getContactSurfaceEstimate() const;

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
  int num_contacts_;
};

} // namespace inekf

#endif // INEKF_SLIP_ESTIMATOR_HPP_