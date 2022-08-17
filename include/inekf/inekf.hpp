/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   inekf.hpp
 *  @author Ross Hartley
 *  @brief  Header file for Invariant EKF 
 *  @date   September 25, 2018
 **/

#ifndef INEKF_INEKF_HPP_
#define INEKF_INEKF_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "Eigen/Core"
#include "Eigen/LU"
#include "unsupported/Eigen/MatrixFunctions"

#include "inekf/inekf_state.hpp"
#include "inekf/noise_params.hpp"
#include "inekf/lie_group.hpp"
#include "inekf/observations.hpp"
#include "inekf/macros.hpp"
#include "inekf/error_type.hpp"
#include "inekf/state_transition_matrix.hpp"
#include "inekf/discrete_noise_matrix.hpp"


namespace inekf {

class InEKF {
public:
/// @name Constructors
/// @{
  /**
   * Default Constructor. Initializes the filter with default state (identity rotation, zero velocity, zero position) and noise parameters. 
   * No contacts, prior landmarks, or magnetic field is set.            gfc = [

    */
  InEKF();
  /**
   * Initialize filter with noise parameters. Initializes th            gfc = [
the default (identity rotation, zero velocity, zero position).
    * @param params: The noise parameters to be assigned.
    */
  InEKF(const NoiseParams& params);
  /**
   * Initialize filter with state. Initializes the noise par            gfc = [
the default.
    * @param state: The state to be assigned.
    */
  InEKF(const InEKFState& state);
  /**
   * Initialize filter with state and noise parameters.
   * @param state: The state to be assigned.
   * @param params: The noise parameters to be assigned.
   */        
  InEKF(const InEKFState& state, const NoiseParams& params);
  /**
   * Initialize filter with state, noise, and error type.
   * @param state: The state to be assigned.
   * @param params: The noise parameters to be assigned.
   * @param error_type: The type of invariant error to be used (affects covariance).
   */       
  InEKF(const InEKFState& state, const NoiseParams& params, const ErrorType error_type);

  INEKF_USE_DEFAULT_DESTTUCTOR(InEKF);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(InEKF);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(InEKF);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(InEKF);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(InEKF);
/// @}

/// @name Getters
/// @{
  /**
   * Gets the current error type.
   */
  ErrorType getErrorType() const;
  /**
   * Gets the current state estimate.
   */
  const InEKFState& getState() const;
  /**
   * Gets the current noise parameters.
   */
  const NoiseParams& getNoiseParams() const;
  /**
   * Gets the filter's current contact states.
   * @return  map of contact ID and bool that indicates if contact is registed
   */
  const std::map<int, bool>& getContacts() const;
  /**
   * Gets the current estimated contact positions.
   * @return  map of contact ID and associated index in the state matrix X
   */
  const std::map<int, int>& getEstimatedContactPositions() const;

  /**
   * Gets the filter's prior landmarks.
   * @return  map of prior landmark ID and position (as a Eigen::Vector3d)
   */
  const mapIntVector3d& getPriorLandmarks() const;
  /**
   * Gets the filter's estimated landmarks.
   * @return  map of landmark ID and associated index in the state matrix X
   */
  const std::map<int, int>& getEstimatedLandmarks() const;
  /**
   * Gets the filter's set magnetic field.
   * @return  magnetic field in world frame
   */
  const Eigen::Vector3d& getMagneticField() const;
/// @}


/// @name Setters
/// @{
  /**
   * Sets the current state estimate
   * @param state: The state estimate to be assigned.
   */
  void setState(const InEKFState& state);
  /**
   * Sets the current noise parameters
   * @param params: The noise parameters to be assigned.
   */
  void setNoiseParams(const NoiseParams& params);
  /**
   * Sets the filter's current contact state.
   * @param contacts: A vector of contact ID and indicator pairs. A true indicator means contact is detected.
   */
  void setContacts(const std::vector<std::pair<int,bool>>& contacts);
  /**
   * Sets the filter's prior landmarks.
   * @param prior_landmarks: A map of prior landmark IDs and associated position in the world frame.
   */
  void setPriorLandmarks(const mapIntVector3d& prior_landmarks);
  /** TODO: Sets magnetic field for untested magnetometer measurement */
  void setMagneticField(const Eigen::Vector3d& true_magnetic_field);
/// @}


/// @name Basic Utilities
/// @{
  /**
   * Resets the filter
   * Initializes state matrix to identity, removes all augmented states, and assigns default noise parameters.
   */
  void clear();
  /**
   * Removes a single landmark from the filter's prior landmark set.
   * @param landmark_id: The ID for the landmark to remove.
   */
  void RemovePriorLandmarks(const int landmark_id);
  /**
   * Removes a set of landmarks from the filter's prior landmark set.
   * @param landmark_ids: A vector of IDs for the landmarks to remove.
   */
  void RemovePriorLandmarks(const std::vector<int>& landmark_ids);
  /**
   * Removes a single landmark from the filter's estimated landmark set.
   * @param landmark_id: The ID for the landmark to remove.
   */
  void RemoveLandmarks(const int landmark_id);
  /**
   * Removes a set of landmarks from the filter's estimated landmark set.
   * @param landmark_ids: A vector of IDs for the landmarks to remove.
   */
  void RemoveLandmarks(const std::vector<int>& landmark_ids);
  /**
   * Keeps a set of landmarks from the filter's estimated landmark set.
   * @param landmark_ids: A vector of IDs for the landmarks to keep.
   */
  void KeepLandmarks(const std::vector<int>& landmark_ids);
/// @}


/// @name Propagation and Correction Methods
/// @{
  /**
   * Propagates the estimated state mean and covariance forward using inertial measurements. 
   * All landmarks positions are assumed to be static.
   * All contacts velocities are assumed to be zero + Gaussian noise.
   * The propagation model currently assumes that the covariance is for the right invariant error.
   * @param imu_w: IMU angular velocity measurement
   * @param imu_a: IMU linear acceleration measurement
   * @param dt: double indicating how long to integrate the inertial measurements for
   */
  void Propagate(const Eigen::Vector3d& imu_w, const Eigen::Vector3d& imu_a, const double dt);
  /**
   * Propagates the estimated state mean and covariance forward using inertial measurements. 
   * All landmarks positions are assumed to be static.
   * All contacts velocities are assumed to be zero + Gaussian noise.
   * The propagation model currently assumes that the covariance is for the right invariant error.
   * @param imu: 6x1 vector containing stacked angular velocity and linear acceleration measurements
   * @param dt: double indicating how long to integrate the inertial measurements for
   */
  void Propagate(const Eigen::Matrix<double,6,1>& imu, const double dt);
  /** 
   * Corrects the state estimate using the measured forward kinematics between the IMU and a set of contact frames.
   * If contact is indicated but not included in the state, the state is augmented to include the estimated contact position.
   * If contact is not indicated but is included in the state, the contact position is marginalized out of the state. 
   * This is a right-invariant measurement model. Example usage can be found in @include kinematics.cpp
   * @param measured_kinematics: the measured kinematics containing the contact id, relative pose measurement in the IMU frame, and covariance
   */
  void CorrectKinematics(const vectorKinematics& measured_kinematics); 
  /** 
   * Corrects the state estimate using the measured position between a set of contact frames and the IMU.
   * If the landmark is not included in the state, the state is augmented to include the estimated landmark position. 
   * This is a right-invariant measurement model.
   * @param measured_landmarks: the measured landmarks containing the contact id, relative position measurement in the IMU frame, and covariance
   */
  void CorrectLandmarks(const vectorLandmarks& measured_landmarks);

  /** TODO: Untested magnetometer measurement*/
  void CorrectMagnetometer(const Eigen::Vector3d& measured_magnetic_field, const Eigen::Matrix3d& covariance);
  /** TODO: Untested GPS measurement*/
  void CorrectPosition(const Eigen::Vector3d& measured_position, const Eigen::Matrix3d& covariance, const Eigen::Vector3d& indices);
  /** TODO: Untested contact position measurement*/
  void CorrectContactPosition(const int id, const Eigen::Vector3d& measured_contact_position, const Eigen::Matrix3d& covariance, const Eigen::Vector3d& indices);
/// @} 

/** @example kinematics.cpp
 * Testing
 */

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ErrorType error_type_ = ErrorType::LeftInvariant; 
  bool estimate_bias_ = true;  
  InEKFState state_;
  NoiseParams noise_params_;
  Eigen::Vector3d g_; // Gravity vector in world frame (z-up)
  std::map<int,bool> contacts_;
  std::map<int,int> estimated_contact_positions_;
  mapIntVector3d prior_landmarks_;
  std::map<int,int> estimated_landmarks_;
  Eigen::Vector3d magnetic_field_;
  Eigen::LDLT<Eigen::MatrixXd> ldlt_;
  StateTransitionMatrix state_transition_matrix_;
  DiscreteNoiseMatrix discrete_noise_matrix_;
  Eigen::MatrixXd P_pred_, X_pred_;
  Eigen::Vector3d phi_;
  Eigen::Matrix3d G0_, G1_, G2_;

  // Corrects state using invariant observation models
  void CorrectRightInvariant(const Observation& obs);
  void CorrectLeftInvariant(const Observation& obs);
  void CorrectRightInvariant(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N);
  void CorrectLeftInvariant(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N);
  // void CorrectFullState(const Observation& obs); // TODO
};

} // namespace inekf 

#endif // INEKF_INEKF_HPP_
