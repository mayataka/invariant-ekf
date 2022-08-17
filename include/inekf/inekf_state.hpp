/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   inekf_state.hpp
 *  @author Ross Hartley
 *  @brief  Header file for InEKFState
 *  @date   September 25, 2018
 **/

#ifndef INEKF_ROBOTSTATE_HPP_
#define INEKF_ROBOTSTATE_HPP_

#include <iostream>

#include "Eigen/Core"

#include "inekf/macros.hpp"


namespace inekf {

enum StateType {WorldCentric, BodyCentric};

class InEKFState {
public:
  InEKFState();
  InEKFState(const Eigen::MatrixXd& X);
  InEKFState(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta);
  InEKFState(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P);

  INEKF_USE_DEFAULT_DESTTUCTOR(InEKFState);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(InEKFState);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(InEKFState);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(InEKFState);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(InEKFState);

  const Eigen::MatrixXd& getX() const;
  const Eigen::VectorXd& getTheta() const;
  const Eigen::MatrixXd& getP() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 3> getRotation() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 1> getVelocity() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 1> getPosition() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 1> getVector(int id) const;
  const Eigen::VectorBlock<const Eigen::VectorXd, 3> getGyroscopeBias() const;
  const Eigen::VectorBlock<const Eigen::VectorXd, 3> getAccelerometerBias() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 3> getRotationCovariance() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 3> getVelocityCovariance() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 3> getPositionCovariance() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 3> getGyroscopeBiasCovariance() const;
  const Eigen::Block<const Eigen::MatrixXd, 3, 3> getAccelerometerBiasCovariance() const;
  int dimX() const;
  int dimTheta() const;
  int dimP() const;
  const StateType getStateType() const;
  const Eigen::MatrixXd getWorldX() const;
  const Eigen::Matrix3d getWorldRotation() const;
  const Eigen::Vector3d getWorldVelocity() const;
  const Eigen::Vector3d getWorldPosition() const;
  const Eigen::MatrixXd getBodyX() const;
  const Eigen::Matrix3d getBodyRotation() const;
  const Eigen::Vector3d getBodyVelocity() const;
  const Eigen::Vector3d getBodyPosition() const;

  void setX(const Eigen::MatrixXd& X);
  void setP(const Eigen::MatrixXd& P);
  void setTheta(const Eigen::VectorXd& Theta);
  void setRotation(const Eigen::Matrix3d& R);
  void setVelocity(const Eigen::Vector3d& v);
  void setPosition(const Eigen::Vector3d& p);
  void setGyroscopeBias(const Eigen::Vector3d& bg);
  void setAccelerometerBias(const Eigen::Vector3d& ba);
  void setRotationCovariance(const Eigen::Matrix3d& cov);
  void setVelocityCovariance(const Eigen::Matrix3d& cov);
  void setPositionCovariance(const Eigen::Matrix3d& cov);
  void setGyroscopeBiasCovariance(const Eigen::Matrix3d& cov);
  void setAccelerometerBiasCovariance(const Eigen::Matrix3d& cov);
  void copyDiagX(const int n, Eigen::MatrixXd& BigX) const;
  void copyDiagXinv(const int n, Eigen::MatrixXd& BigXinv) const;

  Eigen::MatrixXd calcXinv() const;

  friend std::ostream& operator<<(std::ostream& os, const InEKFState& s);  

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  StateType state_type_ = StateType::WorldCentric; 
  Eigen::MatrixXd X_;
  Eigen::VectorXd Theta_;
  Eigen::MatrixXd P_;
};

} // namespace inekf 

#endif // INEKF_ROBOTSTATE_HPP_