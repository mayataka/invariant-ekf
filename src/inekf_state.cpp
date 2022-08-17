/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley <m.ross.hartley@gmail.com>
 * All Rights Reserved 
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   inekf_state.hpp
 *  @author Ross Hartley
 *  @brief  Source file for InEKFState (thread-safe)
 *  @date   September 25, 2018
 **/

#include "legged_state_estimator/inekf/inekf_state.hpp"
#include "legged_state_estimator/inekf/lie_group.hpp"


namespace legged_state_estimator {

// Default constructor
InEKFState::InEKFState() : 
  X_(Eigen::MatrixXd::Identity(5,5)), 
  Theta_(Eigen::MatrixXd::Zero(6,1)), 
  P_(Eigen::MatrixXd::Identity(15,15)) {}


// Initialize with X
InEKFState::InEKFState(const Eigen::MatrixXd& X) : 
    X_(X), Theta_(Eigen::MatrixXd::Zero(6,1)) {
    P_ = Eigen::MatrixXd::Identity(3*this->dimX()+this->dimTheta()-6, 3*this->dimX()+this->dimTheta()-6);
}


// Initialize with X and Theta
InEKFState::InEKFState(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta) : 
    X_(X), Theta_(Theta) {
    P_ = Eigen::MatrixXd::Identity(3*this->dimX()+this->dimTheta()-6, 3*this->dimX()+this->dimTheta()-6);
}


// Initialize with X, Theta and P
InEKFState::InEKFState(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : 
    X_(X), Theta_(Theta), P_(P) {}
// TODO: error checking to make sure dimensions are correct and supported


const Eigen::MatrixXd& InEKFState::getX() const { return X_; }
const Eigen::VectorXd& InEKFState::getTheta() const { return Theta_; }
const Eigen::MatrixXd& InEKFState::getP() const { return P_; }
const Eigen::Block<const Eigen::MatrixXd, 3, 3> InEKFState::getRotation() const { return X_.template block<3,3>(0,0); }
const Eigen::Block<const Eigen::MatrixXd, 3, 1> InEKFState::getVelocity() const { return X_.template block<3,1>(0,3); }
const Eigen::Block<const Eigen::MatrixXd, 3, 1> InEKFState::getPosition() const { return X_.template block<3,1>(0,4); }
const Eigen::Block<const Eigen::MatrixXd, 3, 1> InEKFState::getVector(int index) const { return X_.template block<3,1>(0,index); }

const Eigen::VectorBlock<const Eigen::VectorXd, 3> InEKFState::getGyroscopeBias() const { return Theta_.template head<3>(); }
const Eigen::VectorBlock<const Eigen::VectorXd, 3> InEKFState::getAccelerometerBias() const { return Theta_.template tail<3>(3); }

const Eigen::Block<const Eigen::MatrixXd, 3, 3> InEKFState::getRotationCovariance() const { return P_.template block<3,3>(0,0); }
const Eigen::Block<const Eigen::MatrixXd, 3, 3> InEKFState::getVelocityCovariance() const { return P_.template block<3,3>(3,3); }
const Eigen::Block<const Eigen::MatrixXd, 3, 3> InEKFState::getPositionCovariance() const { return P_.template block<3,3>(6,6); }
const Eigen::Block<const Eigen::MatrixXd, 3, 3> InEKFState::getGyroscopeBiasCovariance() const { return P_.template block<3,3>(P_.rows()-6,P_.rows()-6); }
const Eigen::Block<const Eigen::MatrixXd, 3, 3> InEKFState::getAccelerometerBiasCovariance() const { return P_.template block<3,3>(P_.rows()-3,P_.rows()-3); }

int InEKFState::dimX() const { return X_.cols(); }
int InEKFState::dimTheta() const {return Theta_.rows();}
int InEKFState::dimP() const { return P_.cols(); }


const StateType InEKFState::getStateType() const { return state_type_; }


const Eigen::MatrixXd InEKFState::getWorldX() const {
  if (state_type_ == StateType::WorldCentric) {
    return this->getX();
  } 
  else {
    return this->calcXinv();
  }
}


const Eigen::Matrix3d InEKFState::getWorldRotation() const {
  if (state_type_ == StateType::WorldCentric) {
    return this->getRotation();
  } 
  else {
    return this->getRotation().transpose();
  }
}


const Eigen::Vector3d InEKFState::getWorldVelocity() const {
  if (state_type_ == StateType::WorldCentric) {
    return this->getVelocity();
  } 
  else {
    return -this->getRotation().transpose()*this->getVelocity();
  }
}


const Eigen::Vector3d InEKFState::getWorldPosition() const {
  if (state_type_ == StateType::WorldCentric) {
    return this->getPosition();
  } 
  else {
    return -this->getRotation().transpose()*this->getPosition();
  }
}


const Eigen::MatrixXd InEKFState::getBodyX() const {
  if (state_type_ == StateType::BodyCentric) {
    return this->getX();
  } 
  else {
    return this->calcXinv();
  }
}


const Eigen::Matrix3d InEKFState::getBodyRotation() const {
  if (state_type_ == StateType::BodyCentric) {
    return this->getRotation();
  } 
  else {
    return this->getRotation().transpose();
  }
}


const Eigen::Vector3d InEKFState::getBodyVelocity() const {
  if (state_type_ == StateType::BodyCentric) {
    return this->getVelocity();
  } 
  else {
    return -this->getRotation().transpose()*this->getVelocity();
  }
}


const Eigen::Vector3d InEKFState::getBodyPosition() const {
  if (state_type_ == StateType::BodyCentric) {
    return this->getPosition();
  } 
  else {
    return -this->getRotation().transpose()*this->getPosition();
  }
}


void InEKFState::setX(const Eigen::MatrixXd& X) { X_ = X; }
void InEKFState::setTheta(const Eigen::VectorXd& Theta) { Theta_ = Theta; }
void InEKFState::setP(const Eigen::MatrixXd& P) { P_ = P; }
void InEKFState::setRotation(const Eigen::Matrix3d& R) { X_.block<3,3>(0,0) = R; }
void InEKFState::setVelocity(const Eigen::Vector3d& v) { X_.block<3,1>(0,3) = v; }
void InEKFState::setPosition(const Eigen::Vector3d& p) { X_.block<3,1>(0,4) = p; }

void InEKFState::setGyroscopeBias(const Eigen::Vector3d& bg) { Theta_.head(3) = bg; }
void InEKFState::setAccelerometerBias(const Eigen::Vector3d& ba) { Theta_.tail(3) = ba; }

void InEKFState::setRotationCovariance(const Eigen::Matrix3d& cov) { P_.block<3,3>(0,0) = cov; }
void InEKFState::setVelocityCovariance(const Eigen::Matrix3d& cov) { P_.block<3,3>(3,3) = cov; }
void InEKFState::setPositionCovariance(const Eigen::Matrix3d& cov) { P_.block<3,3>(6,6) = cov; }
void InEKFState::setGyroscopeBiasCovariance(const Eigen::Matrix3d& cov) { P_.block<3,3>(P_.rows()-6,P_.rows()-6) = cov; }
void InEKFState::setAccelerometerBiasCovariance(const Eigen::Matrix3d& cov) { P_.block<3,3>(P_.rows()-3,P_.rows()-3) = cov; }


void InEKFState::copyDiagX(const int n, Eigen::MatrixXd& BigX) const {
  const int dimX = this->dimX();
  for(int i=0; i<n; ++i) {
    const int startIndex = BigX.rows();
    BigX.conservativeResize(startIndex+dimX, startIndex+dimX);
    BigX.block(startIndex,0,dimX,startIndex).setZero();
    BigX.block(0,startIndex,startIndex,dimX).setZero();
    BigX.block(startIndex,startIndex,dimX,dimX) = X_;
  }
  return;
}


void InEKFState::copyDiagXinv(const int n, Eigen::MatrixXd& BigXinv) const {
  const int dimX = this->dimX();
  Eigen::MatrixXd Xinv = this->calcXinv();
  for(int i=0; i<n; ++i) {
    int startIndex = BigXinv.rows();
    BigXinv.conservativeResize(startIndex + dimX, startIndex + dimX);
    BigXinv.block(startIndex,0,dimX,startIndex).setZero();
    BigXinv.block(0,startIndex,startIndex,dimX).setZero();
    BigXinv.block(startIndex,startIndex,dimX,dimX) = Xinv;
  }
  return;
}


Eigen::MatrixXd InEKFState::calcXinv() const {
  const int dimX = this->dimX();
  Eigen::MatrixXd Xinv = Eigen::MatrixXd::Identity(dimX,dimX);
  const auto& RT = X_.block<3,3>(0,0).transpose();
  Xinv.block<3,3>(0,0) = RT;
  for(int i=3; i<dimX; ++i) {
    Xinv.block<3,1>(0,i).noalias() = -RT * X_.block<3,1>(0,i);
  }
  return Xinv;
}


std::ostream& operator<<(std::ostream& os, const InEKFState& s) {  
  os << "--------- InEKF State -------------" << std::endl;
  os << "X:\n" << s.X_ << std::endl << std::endl;
  os << "Theta:\n" << s.Theta_ << std::endl << std::endl;
  // os << "P:\n" << s.P_ << endl;
  os << "-----------------------------------";
  return os;  
} 

} // namespace legged_state_estimator 
