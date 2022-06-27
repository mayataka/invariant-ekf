#include "inekf/slip_estimator.hpp"


namespace inekf {

SlipEstimator::SlipEstimator(const RobotModel& robot_model, 
                             const SlipEstimatorSettings& settings)
  : settings_(settings),
    contact_surface_estimate_(robot_model.numContacts(), Eigen::Matrix3d::Identity()),
    contact_velocity_(robot_model.numContacts(), Eigen::Vector3d::Zero()), 
    contact_velocity_prev_(robot_model.numContacts(), Eigen::Vector3d::Zero()),
    force_velocity_plane_normal_estimate_(robot_model.numContacts(), Eigen::Vector3d::Zero()),
    contact_surface_normal_estimate_(robot_model.numContacts(), (Eigen::Vector3d() << 0., 0., 1.0).finished()),
    contact_force_surface_local_(robot_model.numContacts(), Eigen::Vector3d::Zero()),
    contact_velocity_norm_(robot_model.numContacts(), 0), 
    slip_probability_(robot_model.numContacts(), 0), 
    slip_covariance_(robot_model.numContacts(), 0), 
    friction_coefficient_estimate_(robot_model.numContacts(), 0),
    slip_state_(),
    num_contacts_(robot_model.numContacts()) {
}


SlipEstimator::SlipEstimator() {
}


void SlipEstimator::reset() {
}


void SlipEstimator::update(const RobotModel& robot_model, 
                           const ContactEstimator& contact_estimator) {
  slip_state_.clear();
  // Estimate slip velocity from robot's differential kinematics
  for (int i=0; i<robot_model.numContacts(); ++i) {
    contact_velocity_[i] = robot_model.getContactVelocityWorld(i);
    contact_velocity_norm_[i] = contact_velocity_[i].norm();
  }
  // Slip probability 
  for (int i=0; i<num_contacts_; ++i) {
    slip_probability_[i]
        = 1.0 / (1.0 + std::exp(- settings_.beta1[i] * contact_velocity_[i].template lpNorm<2>()
                                - settings_.beta0[i]));
    if (std::isnan(slip_probability_[i]) || std::isinf(slip_probability_[i])) {
      slip_probability_[i] = 0;
    }
  }
  // Slip covariance
  for (int i=0; i<num_contacts_; ++i) {
    const double contact_velocity_diff = (contact_velocity_[i] - contact_velocity_prev_[i]).template lpNorm<2>();
    slip_covariance_[i] = settings_.slip_velocity_cov_alpha * contact_velocity_diff * contact_velocity_diff;
    contact_velocity_prev_[i] = contact_velocity_[i];
  }
  // Deterministic slip state
  const auto& contact_state = contact_estimator.getContactState();
  slip_state_.clear();
  for (int i=0; i<num_contacts_; ++i) {
    slip_state_.push_back(
        std::pair<int, bool>(i, 
          (slip_probability_[i] >= settings_.slip_prob_threshold) && contact_state[i].second));
  }
  // Estimate the contact surface and normal
  for (int i=0; i<num_contacts_; ++i) {
    force_velocity_plane_normal_estimate_[i] 
        = contact_velocity_[i].cross(contact_estimator.getContactForceEstimate()[i]).normalized();
    contact_surface_normal_estimate_[i] 
        = force_velocity_plane_normal_estimate_[i].cross(contact_velocity_[i]).normalized();
    const double nx = contact_surface_normal_estimate_[i].coeff(0);
    const double ny = contact_surface_normal_estimate_[i].coeff(1);
    const double nz = contact_surface_normal_estimate_[i].coeff(2);
    const double nxny_norm = std::sqrt(nx*nx + ny*ny);
    contact_surface_estimate_[i]
        <<     ny/nxny_norm,   -nx/nxny_norm,         0.,
            nx*nz/nxny_norm, ny*nz/nxny_norm, -nxny_norm,
                         nx,              ny,         nz;
  }
  // Friction coefficient 
  for (int i=0; i<num_contacts_; ++i) {
    contact_force_surface_local_[i].noalias()
        = contact_surface_estimate_[i].transpose() * contact_estimator.getContactForceEstimate()[i];
    const double fn = contact_force_surface_local_[i].coeff(2);
    const double ft = contact_force_surface_local_[i].template head<2>().template lpNorm<2>();
    friction_coefficient_estimate_[i] = ft / fn;
  }
}


void SlipEstimator::setParameters(const SlipEstimatorSettings& settings) {
  settings_ = settings;
}


const std::vector<std::pair<int, bool>>& SlipEstimator::getSlipState() const {
  return slip_state_; 
}


const std::vector<double>& SlipEstimator::getSlipProbability() const {
  return slip_probability_;
}


const std::vector<double>& SlipEstimator::getSlipVelocityCovariance() const {
  return slip_covariance_;  
}


const std::vector<double>& SlipEstimator::getFrictionCoefficientEstimate() const {
 return friction_coefficient_estimate_;
}


const std::vector<Eigen::Vector3d>& SlipEstimator::getContactSurfaceNormalEstimate() const {
  return contact_surface_normal_estimate_;
}


const std::vector<Eigen::Matrix3d>& SlipEstimator::getContactSurfaceEstimate() const {
  return contact_surface_estimate_;
}

} // namespace inekf