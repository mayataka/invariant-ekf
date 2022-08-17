#ifndef INEKF_CONTACT_ESTIMATOR_HPP_
#define INEKF_CONTACT_ESTIMATOR_HPP_

#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/StdVector"

#include "legged_state_estimator/macros.hpp"
#include "legged_state_estimator/robot_model.hpp"


namespace legged_state_estimator {

struct ContactEstimatorSettings {
  std::vector<double> beta0;
  std::vector<double> beta1;
  double contact_force_cov_alpha;
  double contact_prob_threshold;
};


class ContactEstimator {
public:
  ContactEstimator(const RobotModel& robot_model, 
                   const ContactEstimatorSettings& settings);

  ContactEstimator();

  INEKF_USE_DEFAULT_DESTTUCTOR(ContactEstimator);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(ContactEstimator);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(ContactEstimator);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(ContactEstimator);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(ContactEstimator);

  void reset();

  // M(q) a + h (q, v) = S^T u + J^T f

  void update(const RobotModel& robot_model, const Eigen::VectorXd& tauJ);

  void setParameters(const ContactEstimatorSettings& settings);

  const std::vector<std::pair<int, bool>>& getContactState() const;

  const std::vector<double>& getContactProbability() const;

  const std::vector<double>& getContactForceCovariance() const;

  // double getContactForceCovariance() const;

  const std::vector<Eigen::Vector3d>& getContactForceEstimate() const;

  const std::vector<double>& getNormalContactForceEstimate() const;

  const std::vector<Eigen::Vector3d>& getContactSurfaceNormal() const;

  void setContactSurfaceNormal(const std::vector<Eigen::Vector3d>& contact_surface_normal);

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const ContactEstimator& d);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ContactEstimatorSettings settings_;
  std::vector<Eigen::Vector3d> contact_force_estimate_, 
                               contact_force_estimate_prev_, 
                               contact_surface_normal_;
  std::vector<double> normal_contact_force_estimate_, 
                      normal_contact_force_estimate_prev_, 
                      contact_probability_, contact_force_covariance_;
  std::vector<std::pair<int, bool>> contact_state_;
  int num_contacts_;
};

} // namespace legged_state_estimator

#endif // INEKF_CONTACT_ESTIMATOR_HPP_ 