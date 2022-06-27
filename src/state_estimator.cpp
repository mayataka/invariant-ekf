#include "inekf/state_estimator.hpp"


namespace inekf {

StateEstimator::StateEstimator(const StateEstimatorSettings& settings)
  : settings_(settings),
    inekf_(settings.noise_params),
    leg_kinematics_(),
    robot_model_(settings.path_to_urdf, settings.imu_frame, settings.contact_frames),
    contact_estimator_(robot_model_, settings.contact_estimator_settings),
    lpf_gyro_accel_world_(settings.dt, settings.lpf_gyro_accel_cutoff),
    lpf_lin_accel_world_(settings.dt, settings.lpf_lin_accel_cutoff),
    lpf_dqJ_(settings.dt, settings.lpf_dqJ_cutoff, robot_model_.nJ()),
    lpf_ddqJ_(settings.dt, settings.lpf_ddqJ_cutoff, robot_model_.nJ()),
    lpf_tauJ_(settings.dt, settings.lpf_tauJ_cutoff, robot_model_.nJ()),
    imu_gyro_local_(Vector3d::Zero()), 
    imu_gyro_world_(Vector3d::Zero()), 
    imu_gyro_world_prev_(Vector3d::Zero()), 
    imu_gyro_accel_local_(Vector3d::Zero()), 
    imu_gyro_accel_world_(Vector3d::Zero()), 
    imu_lin_accel_local_(Vector3d::Zero()),
    imu_lin_accel_world_(Vector3d::Zero()), 
    quat_(Eigen::Quaterniond::Identity().coeffs()) {
  const double contact_position_cov = settings.contact_position_noise * settings.contact_position_noise;
  const double contact_rotation_cov = settings.contact_rotation_noise * settings.contact_rotation_noise;
  Matrix6d cov_leg = Matrix6d::Zero();
  cov_leg.topLeftCorner<3, 3>() = contact_position_cov * Eigen::Matrix3d::Identity();
  cov_leg.bottomRightCorner<3, 3>() = contact_rotation_cov * Eigen::Matrix3d::Identity();
  for (int i=0; i<settings.contact_frames.size(); ++i) {
    leg_kinematics_.emplace_back(i, Eigen::Matrix4d::Identity(), cov_leg);
  }
}


StateEstimator::StateEstimator() 
  : settings_(),
    inekf_(),
    leg_kinematics_(),
    robot_model_(),
    contact_estimator_(),
    lpf_gyro_accel_world_(),
    lpf_lin_accel_world_(),
    lpf_dqJ_(),
    lpf_ddqJ_(),
    lpf_tauJ_(),
    imu_gyro_local_(Vector3d::Zero()), 
    imu_gyro_world_(Vector3d::Zero()), 
    imu_gyro_world_prev_(Vector3d::Zero()), 
    imu_gyro_accel_local_(Vector3d::Zero()), 
    imu_gyro_accel_world_(Vector3d::Zero()), 
    imu_lin_accel_local_(Vector3d::Zero()),
    imu_lin_accel_world_(Vector3d::Zero()), 
    quat_(Eigen::Quaterniond::Identity().coeffs()) {
}


StateEstimator::~StateEstimator() {}


void StateEstimator::init(const Eigen::Vector3d& base_pos,
                          const Eigen::Vector4d& base_quat,
                          const Eigen::Vector3d& base_lin_vel_world,
                          const Eigen::Vector3d& imu_gyro_bias,
                          const Eigen::Vector3d& imu_lin_accel_bias) {
  RobotState initial_state;
  initial_state.setPosition(base_pos);
  initial_state.setRotation(Eigen::Quaterniond(base_quat).toRotationMatrix());
  initial_state.setVelocity(base_lin_vel_world);
  initial_state.setGyroscopeBias(imu_gyro_bias);
  initial_state.setAccelerometerBias(imu_lin_accel_bias);
  inekf_.setState(initial_state);
}


void StateEstimator::init(const Eigen::Vector3d& base_pos,
                          const Eigen::Vector4d& base_quat,
                          const Eigen::VectorXd& qJ, 
                          const std::vector<double>& ground_height,
                          const Eigen::Vector3d& base_lin_vel_world,
                          const Eigen::Vector3d& imu_gyro_bias,
                          const Eigen::Vector3d& imu_lin_accel_bias) {
  robot_model_.updateKinematics(Eigen::Vector3d::Zero(), base_quat, qJ);
  double base_height = 0;
  for (int i=0; i<robot_model_.numContacts(); ++i) {
    base_height += (robot_model_.getBasePosition().coeff(2)
                      - robot_model_.getContactPosition(i).coeff(2)
                      + ground_height[i]);
  }
  base_height /= robot_model_.numContacts();
  const Eigen::Vector3d base_pos_correct 
      = {base_pos.coeff(0), base_pos.coeff(1), base_height};
  init(base_pos_correct, base_quat, base_lin_vel_world, 
       imu_gyro_bias, imu_lin_accel_bias);
}


void StateEstimator::update(const Eigen::Vector3d& imu_gyro_raw, 
                            const Eigen::Vector3d& imu_lin_accel_raw, 
                            const Eigen::VectorXd& qJ, 
                            const Eigen::VectorXd& dqJ, 
                            const Eigen::VectorXd& tauJ) {
  // Process IMU measurements in InEKF
  inekf_.Propagate(imu_gyro_raw, imu_lin_accel_raw, settings_.dt);
  // Process IMU measurements in LPFs (linear acceleration)
  imu_gyro_local_ = imu_gyro_raw - getIMUGyroBiasEstimate();
  imu_gyro_world_.noalias() = getBaseRotationEstimate() * imu_gyro_local_;
  imu_lin_accel_local_ = imu_lin_accel_raw - getIMULinearAccelerationBiasEstimate();
  imu_lin_accel_world_.noalias() = getBaseRotationEstimate() * imu_lin_accel_local_;
  lpf_lin_accel_world_.update(imu_lin_accel_world_);
  imu_lin_accel_local_.noalias() = getBaseRotationEstimate().transpose() * lpf_lin_accel_world_.getEstimate();
  // Process IMU measurements in LPFs (angular acceleration via a finite difference)
  if (settings_.dynamic_contact_estimation) {
    imu_gyro_accel_world_.noalias() = (imu_gyro_world_ - imu_gyro_world_prev_) / settings_.dt;
    lpf_gyro_accel_world_.update(imu_gyro_accel_world_);
    imu_gyro_accel_local_.noalias() = getBaseRotationEstimate().transpose() * lpf_gyro_accel_world_.getEstimate();
    imu_gyro_world_prev_ = imu_gyro_world_;
  }
  // Process joint measurements in LPFs 
  if (settings_.dynamic_contact_estimation) {
    lpf_ddqJ_.update((dqJ-lpf_dqJ_.getEstimate())/settings_.dt);
  }
  lpf_dqJ_.update(dqJ);
  lpf_tauJ_.update(tauJ);
  // Update contact info
  robot_model_.updateLegKinematics(qJ);
  if (settings_.dynamic_contact_estimation) {
    robot_model_.updateDynamics(getBasePositionEstimate(), getBaseQuaternionEstimate(),
                                getBaseLinearVelocityEstimateLocal(), imu_gyro_raw,
                                imu_lin_accel_local_, imu_gyro_accel_local_,
                                qJ, dqJ, lpf_ddqJ_.getEstimate());
  }
  else {
    robot_model_.updateLegDynamics(qJ, dqJ);
  }
  contact_estimator_.update(robot_model_, lpf_tauJ_.getEstimate());
  inekf_.setContacts(contact_estimator_.getContactState());
  const double contact_force_cov = contact_estimator_.getContactForceCovariance();
  for (int i=0; i<robot_model_.numContacts(); ++i) {
    leg_kinematics_[i].setContactPosition(
        robot_model_.getContactPosition(i)-robot_model_.getBasePosition());
    leg_kinematics_[i].setContactPositionCovariance(
        contact_force_cov*Eigen::Matrix3d::Identity());
  }
  // Process kinematics measurements in InEKF
  inekf_.CorrectKinematics(leg_kinematics_);
  quat_ = Eigen::Quaterniond(getBaseRotationEstimate()).coeffs();
}


void StateEstimator::updateContactInfo(const Eigen::VectorXd& qJ, 
                                       const Eigen::VectorXd& dqJ) {
  robot_model_.updateKinematics(getBasePositionEstimate(), 
                                getBaseQuaternionEstimate(),
                                getBaseLinearVelocityEstimateLocal(),
                                getBaseAngularVelocityEstimateLocal(), qJ, dqJ);
}


const Eigen::Block<const Eigen::MatrixXd, 3, 1> 
StateEstimator::getBasePositionEstimate() const {
  return inekf_.getState().getPosition();
}


const Eigen::Block<const Eigen::MatrixXd, 3, 3> 
StateEstimator::getBaseRotationEstimate() const {
  return inekf_.getState().getRotation();
}


const Eigen::Vector4d& StateEstimator::getBaseQuaternionEstimate() const {
  return quat_;
}


const Eigen::Block<const Eigen::MatrixXd, 3, 1> 
StateEstimator::getBaseLinearVelocityEstimateWorld() const {
  return inekf_.getState().getVelocity();
}


const Eigen::Vector3d StateEstimator::getBaseLinearVelocityEstimateLocal() const {
  return getBaseRotationEstimate().transpose() * getBaseLinearVelocityEstimateWorld();
}


const Eigen::Vector3d& StateEstimator::getBaseAngularVelocityEstimateWorld() const {
  return imu_gyro_world_;
}


const Eigen::Vector3d& StateEstimator::getBaseAngularVelocityEstimateLocal() const {
  return imu_gyro_local_;
}


const Eigen::VectorBlock<const Eigen::VectorXd, 3> 
StateEstimator::getIMUGyroBiasEstimate() const {
  return inekf_.getState().getGyroscopeBias();
}


const Eigen::VectorBlock<const Eigen::VectorXd, 3> 
StateEstimator::getIMULinearAccelerationBiasEstimate() const {
  return inekf_.getState().getAccelerometerBias();
}


const Eigen::VectorXd& StateEstimator::getJointVelocityEstimate() const {
  return lpf_dqJ_.getEstimate();
}


const Eigen::VectorXd& StateEstimator::getJointAccelerationEstimate() const {
  return lpf_ddqJ_.getEstimate();
}


const Eigen::VectorXd& StateEstimator::getJointTorqueEstimate() const {
  return lpf_tauJ_.getEstimate();
}


const std::vector<Eigen::Vector3d>& StateEstimator::getContactForceEstimate() const {
  return contact_estimator_.getContactForceEstimate();
}


const std::vector<double>& StateEstimator::getContactProbability() const {
  return contact_estimator_.getContactProbability();
}


const StateEstimatorSettings& StateEstimator::getSettings() const {
  return settings_;
}

} // namespace inekf
