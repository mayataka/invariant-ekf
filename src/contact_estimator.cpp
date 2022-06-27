#include "inekf/contact_estimator.hpp"


namespace inekf {

ContactEstimator::ContactEstimator(const RobotModel& robot_model,
                                   const ContactEstimatorSettings& settings)
  : settings_(settings),
    contact_force_estimate_(robot_model.numContacts(), Eigen::Vector3d::Zero()),
    contact_force_estimate_prev_(robot_model.numContacts(), Eigen::Vector3d::Zero()),
    contact_surface_normal_(robot_model.numContacts(), Eigen::Vector3d::Zero()),
    normal_contact_force_estimate_(robot_model.numContacts(), 0),
    normal_contact_force_estimate_prev_(robot_model.numContacts(), 0),
    contact_probability_(robot_model.numContacts(), 0),
    contact_force_covariance_(robot_model.numContacts(), 0),
    contact_state_(),
    num_contacts_(robot_model.numContacts()) {
  for (auto& e : contact_surface_normal_) {
    e << 0, 0, 1;
  }
  contact_state_.clear();
  for (int i=0; i<num_contacts_; ++i) {
    contact_state_.push_back(std::pair<int, bool>(i, false));
  }
}


ContactEstimator::ContactEstimator() 
  : settings_(),
    contact_force_estimate_(),
    contact_force_estimate_prev_(),
    contact_surface_normal_(),
    normal_contact_force_estimate_(),
    normal_contact_force_estimate_prev_(),
    contact_probability_(),
    contact_force_covariance_(),
    contact_state_(),
    num_contacts_(0) {
}


void ContactEstimator::reset() {
}


void ContactEstimator::update(const RobotModel& robot_model, 
                              const Eigen::VectorXd& tauJ) {
  // Estimate contact force from robot dynamics via logistic regression
  for (int i=0; i<num_contacts_; ++i) {
    contact_force_estimate_[i].noalias() 
        = - robot_model.getJointContactJacobian(i).template block<3, 3>(0, i*3).transpose().inverse() 
            * (tauJ.template segment<3>(3*i)-robot_model.getJointInverseDynamics().template segment<3>(3*i));
    normal_contact_force_estimate_[i] = contact_force_estimate_[i].dot(contact_surface_normal_[i]);
  }
  // Contact probability 
  for (int i=0; i<num_contacts_; ++i) {
    contact_probability_[i]
        = 1.0 / (1.0 + std::exp(- settings_.beta1[i] * normal_contact_force_estimate_[i]
                                - settings_.beta0[i]));
    if (std::isnan(contact_probability_[i]) || std::isinf(contact_probability_[i])) {
      contact_probability_[i] = 0;
    }
  }
  // Contact covariance
  for (int i=0; i<num_contacts_; ++i) {
    const double df = normal_contact_force_estimate_[i] - normal_contact_force_estimate_prev_[i];
    contact_force_covariance_[i] = settings_.contact_force_cov_alpha * df * df;
    contact_force_estimate_prev_[i] = contact_force_estimate_[i];
    normal_contact_force_estimate_prev_[i] = normal_contact_force_estimate_[i];
  }
  // Deterministic contact state
  contact_state_.clear();
  for (int i=0; i<num_contacts_; ++i) {
    contact_state_.push_back(
        std::pair<int, bool>(i, (contact_probability_[i] >= settings_.contact_prob_threshold)));
  }
}


void ContactEstimator::setParameters(const ContactEstimatorSettings& settings) {
  settings_ = settings;
}


const std::vector<std::pair<int, bool>>& ContactEstimator::getContactState() const {
  return contact_state_;
}


const std::vector<double>& ContactEstimator::getContactProbability() const {
  return contact_probability_;
}


const std::vector<double>& ContactEstimator::getContactForceCovariance() const {
  return contact_force_covariance_;
}


// double ContactEstimator::getContactForceCovariance() const {
//   int num_active_contacts = 0;
//   for (const auto e : getContactState()) {
//     if (e.second) {
//       ++num_active_contacts;
//     }
//   }
//   double contact_force_cov = 0;
//   for (int i=0; i<num_contacts_; ++i) {
//     if (getContactState()[i].second) {
//       contact_force_cov += contact_force_covariance_[i];
//     }
//   }
//   if (num_active_contacts > 0) {
//     contact_force_cov *= (1.0/num_active_contacts);
//   }
//   else {
//     contact_force_cov = 0;
//   }
//   return contact_force_cov;
// }


const std::vector<Eigen::Vector3d>& ContactEstimator::getContactForceEstimate() const {
  return contact_force_estimate_;
}


const std::vector<double>& ContactEstimator::getNormalContactForceEstimate() const {
  return normal_contact_force_estimate_;
}


const std::vector<Eigen::Vector3d>& ContactEstimator::getContactSurfaceNormal() const {
  return contact_surface_normal_;
}


void ContactEstimator::setContactSurfaceNormal(
    const std::vector<Eigen::Vector3d>& contact_surface_normal) {
  assert(contact_surface_normal_.size() == contact_surface_normal.size());
  contact_surface_normal_ = contact_surface_normal;
}


void ContactEstimator::disp(std::ostream& os) const {
  os << "Contact estimation:" << std::endl;
  os << "  contact state: [";
  for (int i=0; i<num_contacts_-1; ++i) {
    os << std::boolalpha << getContactState()[i].second << ", ";
  }
  os << std::boolalpha << getContactState()[num_contacts_-1].second << "]" << std::endl;
  /////////////////////////////////////////
  os << "  contact probability: [";
  for (int i=0; i<num_contacts_-1; ++i) {
    os << getContactProbability()[i] << ", ";
  }
  os << getContactProbability()[num_contacts_-1] << "]" << std::endl;
  /////////////////////////////////////////
  os << "  contact force covariance: [";
  for (int i=0; i<num_contacts_-1; ++i) {
    os << getContactForceCovariance()[i] << ", ";
  }
  os << getContactForceCovariance()[num_contacts_-1] << "]" << std::endl;
  /////////////////////////////////////////
  os << "  contact force estimate: [";
  for (int i=0; i<num_contacts_-1; ++i) {
    os << "[" << getContactForceEstimate()[i].transpose() << "], ";
  }
  os << "[" << getContactForceEstimate()[num_contacts_-1].transpose() << "]]" << std::endl;
  /////////////////////////////////////////
  os << "  normal contact force estimate: [";
  for (int i=0; i<num_contacts_-1; ++i) {
    os << getNormalContactForceEstimate()[i] << ", ";
  }
  os << getNormalContactForceEstimate()[num_contacts_-1] << "]" << std::endl;
  /////////////////////////////////////////
  os << "  assumed contact surface normal: [";
  for (int i=0; i<num_contacts_-1; ++i) {
    os << "[" << getContactSurfaceNormal()[i].transpose() << "], ";
  }
  os << "[" << getContactSurfaceNormal()[num_contacts_-1].transpose() << "]]" << std::flush;
}


std::ostream& operator<<(std::ostream& os, const ContactEstimator& e) {
  e.disp(os);
  return os;
}

} // namespace inekf
