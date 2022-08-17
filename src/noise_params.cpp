/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley <m.ross.hartley@gmail.com>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   NoiseParams.cpp
 *  @author Ross Hartley
 *  @brief  Source file for Invariant EKF noise parameter class
 *  @date   September 25, 2018
 **/

#include "legged_state_estimator/inekf/noise_params.hpp"


namespace legged_state_estimator {

// ------------ NoiseParams -------------
// Default Constructor
NoiseParams::NoiseParams() {
  setGyroscopeNoise(0.01);
  setAccelerometerNoise(0.1);
  setGyroscopeBiasNoise(0.00001);
  setAccelerometerBiasNoise(0.0001);
  setContactNoise(0.1);
}

void NoiseParams::setGyroscopeNoise(const double stddev) { Qg_ = stddev*stddev*Eigen::Matrix3d::Identity(); }
void NoiseParams::setGyroscopeNoise(const Eigen::Vector3d& stddev) { Qg_ << stddev(0)*stddev(0),0,0, 0,stddev(1)*stddev(1),0, 0,0,stddev(2)*stddev(2); }
void NoiseParams::setGyroscopeNoise(const Eigen::Matrix3d& cov) { Qg_ = cov; }

void NoiseParams::setAccelerometerNoise(const double stddev) { Qa_ = stddev*stddev*Eigen::Matrix3d::Identity(); }
void NoiseParams::setAccelerometerNoise(const Eigen::Vector3d& stddev) { Qa_ << stddev(0)*stddev(0),0,0, 0,stddev(1)*stddev(1),0, 0,0,stddev(2)*stddev(2); }
void NoiseParams::setAccelerometerNoise(const Eigen::Matrix3d& cov) { Qa_ = cov; } 

void NoiseParams::setGyroscopeBiasNoise(const double stddev) { Qbg_ = stddev*stddev*Eigen::Matrix3d::Identity(); }
void NoiseParams::setGyroscopeBiasNoise(const Eigen::Vector3d& stddev) { Qbg_ << stddev(0)*stddev(0),0,0, 0,stddev(1)*stddev(1),0, 0,0,stddev(2)*stddev(2); }
void NoiseParams::setGyroscopeBiasNoise(const Eigen::Matrix3d& cov) { Qbg_ = cov; }

void NoiseParams::setAccelerometerBiasNoise(const double stddev) { Qba_ = stddev*stddev*Eigen::Matrix3d::Identity(); }
void NoiseParams::setAccelerometerBiasNoise(const Eigen::Vector3d& stddev) { Qba_ << stddev(0)*stddev(0),0,0, 0,stddev(1)*stddev(1),0, 0,0,stddev(2)*stddev(2); }
void NoiseParams::setAccelerometerBiasNoise(const Eigen::Matrix3d& cov) { Qba_ = cov; }

void NoiseParams::setContactNoise(const double stddev) { Qc_ = stddev*stddev*Eigen::Matrix3d::Identity(); }
void NoiseParams::setContactNoise(const Eigen::Vector3d& stddev) { Qc_ << stddev(0)*stddev(0),0,0, 0,stddev(1)*stddev(1),0, 0,0,stddev(2)*stddev(2); }
void NoiseParams::setContactNoise(const Eigen::Matrix3d& cov) { Qc_ = cov; }

const Eigen::Matrix3d& NoiseParams::getGyroscopeCov() const { return Qg_; }
const Eigen::Matrix3d& NoiseParams::getAccelerometerCov() const { return Qa_; }
const Eigen::Matrix3d& NoiseParams::getGyroscopeBiasCov() const { return Qbg_; }
const Eigen::Matrix3d& NoiseParams::getAccelerometerBiasCov() const { return Qba_; }
const Eigen::Matrix3d& NoiseParams::getContactCov() const { return Qc_; }

std::ostream& operator<<(std::ostream& os, const NoiseParams& p) {
  os << "--------- Noise Params -------------" << std::endl;
  os << "Gyroscope Covariance:\n" << p.Qg_ << std::endl;
  os << "Accelerometer Covariance:\n" << p.Qa_ << std::endl;
  os << "Gyroscope Bias Covariance:\n" << p.Qbg_ << std::endl;
  os << "Accelerometer Bias Covariance:\n" << p.Qba_ << std::endl;
  os << "Contact Covariance:\n" << p.Qc_ << std::endl;
  os << "-----------------------------------" << std::endl;
  return os;
}

} // end inekf namespace
