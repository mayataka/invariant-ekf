/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley <m.ross.hartley@gmail.com>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   InEKF.cpp
 *  @author Ross Hartley
 *  @brief  Source file for Invariant EKF 
 *  @date   September 25, 2018
 **/

#include "inekf/inekf.hpp"


namespace inekf {

using namespace std;

void removeRowAndColumn(Eigen::MatrixXd& M, int index);

// Default constructor
InEKF::InEKF() 
  : g_((Eigen::VectorXd(3) << 0,0,-9.81).finished()), 
    magnetic_field_((Eigen::VectorXd(3) << 0,0,0).finished()),
    state_transition_matrix_(),
    discrete_noise_matrix_() {}

// Constructor with noise params
InEKF::InEKF(const NoiseParams& params) 
  : g_((Eigen::VectorXd(3) << 0,0,-9.81).finished()), 
    magnetic_field_((Eigen::VectorXd(3) << std::cos(1.2049),0,std::sin(1.2049)).finished()), 
    noise_params_(params),
    state_transition_matrix_(),
    discrete_noise_matrix_() {}

// Constructor with initial state
InEKF::InEKF(const RobotState& state) 
  : g_((Eigen::VectorXd(3) << 0,0,-9.81).finished()), 
    magnetic_field_((Eigen::VectorXd(3) << std::cos(1.2049),0,std::sin(1.2049)).finished()), 
    state_(state),
    state_transition_matrix_(state),
    discrete_noise_matrix_(state) {
  const int dimX = state.dimX();
  const int dimTheta = state.dimTheta();
  const int dimP = state.dimP();    
  P_pred_.setZero(dimP, dimP);
  X_pred_.setZero(dimX, dimX);
}

// Constructor with initial state and noise params
InEKF::InEKF(const RobotState& state, const NoiseParams& params) 
  : g_((Eigen::VectorXd(3) << 0,0,-9.81).finished()), 
    magnetic_field_((Eigen::VectorXd(3) << std::cos(1.2049),0,std::sin(1.2049)).finished()), 
    state_(state), 
    noise_params_(params),
    state_transition_matrix_(state),
    discrete_noise_matrix_(state) {
  const int dimX = state.dimX();
  const int dimTheta = state.dimTheta();
  const int dimP = state.dimP();    
  P_pred_.setZero(dimP, dimP);
  X_pred_.setZero(dimX, dimX);
}

// Constructor with initial state, noise params, and error type
InEKF::InEKF(const RobotState& state, const NoiseParams& params, const ErrorType error_type) 
  : g_((Eigen::VectorXd(3) << 0,0,-9.81).finished()), 
    magnetic_field_((Eigen::VectorXd(3) << std::cos(1.2049),0,std::sin(1.2049)).finished()), 
    state_(state), 
    noise_params_(params), 
    error_type_(error_type),
    state_transition_matrix_(state, error_type),
    discrete_noise_matrix_(state, error_type) {
  const int dimX = state.dimX();
  const int dimTheta = state.dimTheta();
  const int dimP = state.dimP();    
  P_pred_.setZero(dimP, dimP);
  X_pred_.setZero(dimX, dimX);
}

// Clear all data in the filter
void InEKF::clear() {
  state_ = RobotState();
  noise_params_ = NoiseParams();
  prior_landmarks_.clear();
  estimated_landmarks_.clear();
  contacts_.clear();
  estimated_contact_positions_.clear();
}

// Returns the robot's current error type
ErrorType InEKF::getErrorType() const { return error_type_; }

// Return robot's current state
const RobotState& InEKF::getState() const { return state_; }

// Sets the robot's current state
void InEKF::setState(const RobotState& state) { state_ = state; }

// Return noise params
const NoiseParams& InEKF::getNoiseParams() const { return noise_params_; }

// Sets the filter's noise parameters
void InEKF::setNoiseParams(const NoiseParams& params) { noise_params_ = params; }

// Return filter's prior (static) landmarks
const mapIntVector3d& InEKF::getPriorLandmarks() const { return prior_landmarks_; }

// Set the filter's prior (static) landmarks
void InEKF::setPriorLandmarks(const mapIntVector3d& prior_landmarks) { prior_landmarks_ = prior_landmarks; }

// Return filter's estimated landmarks
const std::map<int,int>& InEKF::getEstimatedLandmarks() const { return estimated_landmarks_; }

// Return filter's estimated landmarks
const std::map<int,int>& InEKF::getEstimatedContactPositions() const { return estimated_contact_positions_; }

// Set the filter's contact state
void InEKF::setContacts(const vector<std::pair<int,bool>>& contacts) {
  // Insert new measured contact states
  for (const auto& e : contacts) {
    std::pair<map<int,bool>::iterator,bool> ret = contacts_.insert(e);
    // If contact is already in the map, replace with new value
    if (ret.second==false) {
      ret.first->second = e.second;
    }
  }
  return;
}

// Return the filter's contact state
const std::map<int,bool>& InEKF::getContacts() const { return contacts_; }

// Set the true magnetic field
void InEKF::setMagneticField(const Eigen::Vector3d& true_magnetic_field) { magnetic_field_ = true_magnetic_field; }

// Get the true magnetic field
const Eigen::Vector3d& InEKF::getMagneticField() const { return magnetic_field_; }


// InEKF Propagation - Inertial Data
void InEKF::Propagate(const Eigen::Vector3d& imu_w, const Eigen::Vector3d& imu_a, double dt) {
  // Bias corrected IMU measurements
  const Eigen::Vector3d w = imu_w - state_.getGyroscopeBias();    // Angular Velocity
  const Eigen::Vector3d a = imu_a - state_.getAccelerometerBias(); // Linear Acceleration

  // Get current state estimate and dimensions
  const auto& X = state_.getX();
  const Eigen::MatrixXd Xinv = state_.calcXinv();
  const auto& P = state_.getP();
  int dimX = state_.dimX();
  int dimP = state_.dimP();
  int dimTheta = state_.dimTheta();
  if (P_pred_.cols() != dimP || P_pred_.rows() != dimP) {
    P_pred_.resize(dimP, dimP);
  }
  if (X_pred_.cols() != dimX || X_pred_.rows() != dimX) {
    X_pred_.resize(dimX, dimX);
  }

  //  ------------ Propagate Covariance --------------- //
  state_transition_matrix_.compute(state_, w, a, dt);
  const Eigen::MatrixXd& Phi = state_transition_matrix_.Phi();
  discrete_noise_matrix_.compute(state_, noise_params_, estimated_contact_positions_, Phi, dt);
  const Eigen::MatrixXd& Qd = discrete_noise_matrix_.Qd();
  P_pred_.noalias() = Phi * P * Phi.transpose() + Qd;

  // If we don't want to estimate bias, remove correlation
  if (!estimate_bias_) {
    P_pred_.block(0,dimP-dimTheta,dimP-dimTheta,dimTheta).setZero();
    P_pred_.block(dimP-dimTheta,0,dimTheta,dimP-dimTheta).setZero();
    P_pred_.block(dimP-dimTheta,dimP-dimTheta,dimTheta,dimTheta).setIdentity();
  }    

  //  ------------ Propagate Mean --------------- // 
  const auto& R = state_.getRotation();
  const auto& v = state_.getVelocity();
  const auto& p = state_.getPosition();
  phi_ = w*dt;
  Gamma_SO3(phi_, G0_, 0); // Computation can be sped up by computing G0,G1,G2 all at once
  Gamma_SO3(phi_, G1_, 1);
  Gamma_SO3(phi_, G2_, 2);

  X_pred_ = X;
  if (state_.getStateType() == StateType::WorldCentric) {
    // Propagate world-centric state estimate
    X_pred_.template block<3,3>(0,0).noalias() = R * G0_;
    X_pred_.template block<3,1>(0,3).noalias() = v + (R*G1_*a + g_)*dt;
    X_pred_.template block<3,1>(0,4).noalias() = p + v*dt + (R*G2_*a + 0.5*g_)*dt*dt;
  } else {
    // Propagate body-centric state estimate
    const auto& G0t = G0_.transpose();
    X_pred_.template block<3,3>(0,0).noalias() = G0t * R;
    X_pred_.template block<3,1>(0,3).noalias() = G0t * (v - (G1_*a + R*g_)*dt);
    X_pred_.template block<3,1>(0,4).noalias() = G0t * (p + v*dt - (G2_*a + 0.5*R*g_)*dt*dt);
    for (int i=5; i<dimX; ++i) {
      X_pred_.template block<3,1>(0,i).noalias() = G0t * X.block<3,1>(0,i);
    }
  } 
  //  ------------ Update State --------------- // 
  state_.setX(X_pred_);
  state_.setP(P_pred_);      
}


void InEKF::Propagate(const Eigen::Matrix<double,6,1>& imu, double dt) {
  Propagate(imu.template head<3>(), imu.template tail<3>(), dt);
}


// Correct State: Right-Invariant Observation
void InEKF::CorrectRightInvariant(const Eigen::MatrixXd& Z, 
                                  const Eigen::MatrixXd& H, 
                                  const Eigen::MatrixXd& N) {
  // Get current state estimate
  const auto& X = state_.getX();
  Eigen::VectorXd Theta = state_.getTheta();
  Eigen::MatrixXd P = state_.getP();
  const int dimX = state_.dimX();
  const int dimTheta = state_.dimTheta();
  const int dimP = state_.dimP();

  // Remove bias
  Theta = Eigen::Matrix<double,6,1>::Zero();
  P.block<6,6>(dimP-dimTheta,dimP-dimTheta) = 0.0001*Eigen::Matrix<double,6,6>::Identity();
  P.block(0,dimP-dimTheta,dimP-dimTheta,dimTheta).setZero();
  P.block(dimP-dimTheta,0,dimTheta,dimP-dimTheta).setZero();
  // std::cout << "P:\n" << P << std::endl;
  // std::cout << state_ << std::endl;

  // Map from left invariant to right invariant error temporarily
  if (error_type_==ErrorType::LeftInvariant) {
    Eigen::MatrixXd Adj = Eigen::MatrixXd::Identity(dimP,dimP);
    Adj.block(0,0,dimP-dimTheta,dimP-dimTheta) = Adjoint_SEK3(X); 
    P.noalias() = (Adj * P * Adj.transpose()).eval(); 
  }

  // Compute Kalman Gain
  const Eigen::MatrixXd PHT = P * H.transpose();
  const Eigen::MatrixXd S = H * PHT + N;
  Eigen::MatrixXd Sinv;
  if (S.rows() <= 3) {
    Sinv = S.inverse();
  }
  else {
    ldlt_.compute(S);
    const int dimS = S.rows();
    Sinv = ldlt_.solve(Eigen::MatrixXd::Identity(dimS, dimS));
  }
  const Eigen::MatrixXd K = PHT * Sinv;

  // Compute state correction vector
  const Eigen::VectorXd delta = K*Z;
  const Eigen::MatrixXd dX = Exp_SEK3(delta.segment(0,delta.rows()-dimTheta));
  const Eigen::VectorXd dTheta = delta.segment(delta.rows()-dimTheta, dimTheta);

  // Update state
  const Eigen::MatrixXd X_new = dX*X; // Right-Invariant Update
  const Eigen::VectorXd Theta_new = Theta + dTheta;

  // Set new state  
  state_.setX(X_new); 
  state_.setTheta(Theta_new);

  // Update Covariance
  const Eigen::MatrixXd IKH = Eigen::MatrixXd::Identity(dimP,dimP) - K*H;
  Eigen::MatrixXd P_new = IKH * P * IKH.transpose() + K*N*K.transpose(); // Joseph update form

  // Map from right invariant back to left invariant error
  if (error_type_==ErrorType::LeftInvariant) {
    Eigen::MatrixXd AdjInv = Eigen::MatrixXd::Identity(dimP,dimP);
    AdjInv.block(0,0,dimP-dimTheta,dimP-dimTheta) = Adjoint_SEK3(state_.calcXinv()); 
    P_new = (AdjInv * P_new * AdjInv.transpose()).eval();
  }
  // Set new covariance
  state_.setP(P_new); 
}   


// Correct State: Left-Invariant Observation
void InEKF::CorrectLeftInvariant(const Eigen::MatrixXd& Z, 
                                 const Eigen::MatrixXd& H, 
                                 const Eigen::MatrixXd& N) {
  // Get current state estimate
  const auto& X = state_.getX();
  const auto& Theta = state_.getTheta();
  Eigen::MatrixXd P = state_.getP();
  int dimX = state_.dimX();
  int dimTheta = state_.dimTheta();
  int dimP = state_.dimP();

  // Map from right invariant to left invariant error temporarily
  if (error_type_==ErrorType::RightInvariant) {
    Eigen::MatrixXd AdjInv = Eigen::MatrixXd::Identity(dimP,dimP);
    AdjInv.block(0,0,dimP-dimTheta,dimP-dimTheta) = Adjoint_SEK3(state_.calcXinv()); 
    P = (AdjInv * P * AdjInv.transpose()).eval();
  }

  // Compute Kalman Gain
  const Eigen::MatrixXd PHT = P * H.transpose();
  const Eigen::MatrixXd S = H * PHT + N;
  Eigen::MatrixXd Sinv;
  if (S.rows() <= 3) {
    Sinv = S.inverse();
  }
  else {
    ldlt_.compute(S);
    const int dimS = S.rows();
    Sinv = ldlt_.solve(Eigen::MatrixXd::Identity(dimS, dimS));
  }
  const Eigen::MatrixXd K = PHT * Sinv;

  // Compute state correction vector
  const Eigen::VectorXd delta = K*Z;
  const Eigen::MatrixXd dX = Exp_SEK3(delta.segment(0,delta.rows()-dimTheta));
  const Eigen::VectorXd dTheta = delta.segment(delta.rows()-dimTheta, dimTheta);

  // Update state
  const Eigen::MatrixXd X_new = X*dX; // Left-Invariant Update
  const Eigen::VectorXd Theta_new = Theta + dTheta;

  // Set new state
  state_.setX(X_new); 
  state_.setTheta(Theta_new);

  // Update Covariance
  const Eigen::MatrixXd IKH = Eigen::MatrixXd::Identity(dimP,dimP) - K*H;
  Eigen::MatrixXd P_new = IKH * P * IKH.transpose() + K*N*K.transpose(); // Joseph update form

  // Map from left invariant back to right invariant error
  if (error_type_==ErrorType::RightInvariant) {
    Eigen::MatrixXd Adj = Eigen::MatrixXd::Identity(dimP,dimP);
    Adj.block(0,0,dimP-dimTheta,dimP-dimTheta) = Adjoint_SEK3(X_new); 
    P_new = (Adj * P_new * Adj.transpose()).eval(); 
  }

  // Set new covariance
  state_.setP(P_new); 
}   

// Correct state using kinematics measured between imu and contact point
void InEKF::CorrectKinematics(const vectorKinematics& measured_kinematics) {
  Eigen::VectorXd Z, Y, b;
  Eigen::MatrixXd H, N, PI;

  vector<pair<int,int> > remove_contacts;
  vectorKinematics new_contacts;
  vector<int> used_contact_ids;

  for (vectorKinematicsIterator it=measured_kinematics.begin(); it!=measured_kinematics.end(); ++it) {
    // Detect and skip if an ID is not unique (this would cause singularity issues in InEKF::Correct)
    if (find(used_contact_ids.begin(), used_contact_ids.end(), it->id) != used_contact_ids.end()) { 
      cout << "Duplicate contact ID detected! Skipping measurement.\n";
      continue; 
    } 
    else { 
      used_contact_ids.push_back(it->id); 
    }

    // Find contact indicator for the kinematics measurement
    map<int,bool>::iterator it_contact = contacts_.find(it->id);
    if (it_contact == contacts_.end()) { continue; } // Skip if contact state is unknown
    bool contact_indicated = it_contact->second;

    // See if we can find id estimated_contact_positions
    map<int,int>::iterator it_estimated = estimated_contact_positions_.find(it->id);
    bool found = it_estimated!=estimated_contact_positions_.end();

    if (!contact_indicated && found) {
      // If contact is not indicated and id is found in estimated_contacts_, then remove state
      remove_contacts.push_back(*it_estimated); // Add id to remove list
    } 
    else if (contact_indicated && !found) {
      // If contact is indicated and id is not found i n estimated_contacts_, then augment state
      new_contacts.push_back(*it); // Add to augment list

    } 
    else if (contact_indicated && found) {
      // If contact is indicated and id is found in estimated_contacts_, then correct using kinematics
      const int dimX = state_.dimX();
      const int dimTheta = state_.dimTheta();
      const int dimP = state_.dimP();
      int startIndex;
      // Fill out H
      startIndex = H.rows();
      H.conservativeResize(startIndex+3, dimP);
      H.block(startIndex,0,3,dimP).setZero();
      if (state_.getStateType() == StateType::WorldCentric) {
        H.template block<3,3>(startIndex,6) = -Eigen::Matrix3d::Identity(); // -I
        H.template block<3,3>(startIndex,3*it_estimated->second-dimTheta) = Eigen::Matrix3d::Identity(); // I
      } 
      else {
        H.template block<3,3>(startIndex,6) = Eigen::Matrix3d::Identity(); // I
        H.template block<3,3>(startIndex,3*it_estimated->second-dimTheta) = -Eigen::Matrix3d::Identity(); // -I
      }
      // Fill out N
      startIndex = N.rows();
      N.conservativeResize(startIndex+3, startIndex+3);
      N.block(startIndex,0,3,startIndex).setZero();
      N.block(0,startIndex,startIndex,3).setZero();
      N.template block<3,3>(startIndex,startIndex).noalias() = state_.getWorldRotation() * it->covariance.block<3,3>(3,3) 
                                                                                          * state_.getWorldRotation().transpose();
      // Fill out Z
      startIndex = Z.rows();
      Z.conservativeResize(startIndex+3, Eigen::NoChange);
      const auto& R = state_.getRotation();
      const auto& p = state_.getPosition();
      const auto& d = state_.getVector(it_estimated->second);  
      if (state_.getStateType() == StateType::WorldCentric) {
        Z.template segment<3>(startIndex).noalias() = R * it->pose.block<3,1>(0,3) - (d - p); 
      } 
      else {
        Z.template segment<3>(startIndex).noalias() = R.transpose() * (it->pose.block<3,1>(0,3) - (p - d)); 
      }
    } 
    else {
      // If contact is not indicated and id is found in estimated_contacts_, then skip
      continue;
    }
  }

  // Correct state using stacked observation
  if (Z.rows()>0) {
    if (state_.getStateType() == StateType::WorldCentric) {
      this->CorrectRightInvariant(Z,H,N);
      // this->CorrectRightInvariant(obs);
    } 
    else {
      // this->CorrectLeftInvariant(obs);
      this->CorrectLeftInvariant(Z,H,N);
    }
  }

  // Remove contacts from state
  if (remove_contacts.size() > 0) {
    Eigen::MatrixXd X_rem = state_.getX(); 
    Eigen::MatrixXd P_rem = state_.getP();
    for (vector<pair<int,int> >::iterator it=remove_contacts.begin(); it!=remove_contacts.end(); ++it) {
      // Remove row and column from X
      removeRowAndColumn(X_rem, it->second);
      // Remove 3 rows and columns from P
      int startIndex = 3 + 3*(it->second-3);
      removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
      removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
      removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
      // Update all indices for estimated_landmarks and estimated_contact_positions
      for (map<int,int>::iterator it2=estimated_landmarks_.begin(); it2!=estimated_landmarks_.end(); ++it2) {
        if (it2->second > it->second) it2->second -= 1;
      }
      for (map<int,int>::iterator it2=estimated_contact_positions_.begin(); it2!=estimated_contact_positions_.end(); ++it2) {
        if (it2->second > it->second) it2->second -= 1;
      }
      // We also need to update the indices of remove_contacts in the case where multiple contacts are being removed at once
      for (vector<pair<int,int> >::iterator it2=it; it2!=remove_contacts.end(); ++it2) {
        if (it2->second > it->second) it2->second -= 1;
      }
      // Remove from list of estimated contact positions 
      estimated_contact_positions_.erase(it->first);
    }
    // Update state and covariance
    state_.setX(X_rem);
    state_.setP(P_rem);
  }


  // Augment state with newly detected contacts
  if (new_contacts.size() > 0) {
    Eigen::MatrixXd X_aug = state_.getX(); 
    Eigen::MatrixXd P_aug = state_.getP();
    for (vectorKinematicsIterator it=new_contacts.begin(); it!=new_contacts.end(); ++it) {
      // Initialize new landmark mean
      int startIndex = X_aug.rows();
      X_aug.conservativeResize(startIndex+1, startIndex+1);
      X_aug.block(startIndex,0,1,startIndex).setZero();
      X_aug.block(0,startIndex,startIndex,1).setZero();
      X_aug(startIndex, startIndex) = 1;
      if (state_.getStateType() == StateType::WorldCentric) {
        X_aug.block(0,startIndex,3,1).noalias() = state_.getPosition() + state_.getRotation() * it->pose.block<3,1>(0,3);
      } 
      else {
        X_aug.block(0,startIndex,3,1).noalias() = state_.getPosition() - it->pose.block<3,1>(0,3);
      }

      // Initialize new landmark covariance - TODO:speed up
      Eigen::MatrixXd F = Eigen::MatrixXd::Zero(state_.dimP()+3,state_.dimP()); 
      F.block(0,0,state_.dimP()-state_.dimTheta(),state_.dimP()-state_.dimTheta()).setIdentity(); // for old X
      F.block(state_.dimP()-state_.dimTheta()+3,state_.dimP()-state_.dimTheta(),state_.dimTheta(),state_.dimTheta()).setIdentity(); // for theta
      Eigen::MatrixXd G = Eigen::MatrixXd::Zero(F.rows(),3);
      // Blocks for new contact
      if ((state_.getStateType() == StateType::WorldCentric && error_type_ == ErrorType::RightInvariant) || 
          (state_.getStateType() == StateType::BodyCentric && error_type_ == ErrorType::LeftInvariant)) {
        F.block(state_.dimP()-state_.dimTheta(),6,3,3) = Eigen::Matrix3d::Identity(); 
        G.block(G.rows()-state_.dimTheta()-3,0,3,3) = state_.getWorldRotation();
      } 
      else {
        F.block(state_.dimP()-state_.dimTheta(),6,3,3) = Eigen::Matrix3d::Identity(); 
        F.block(state_.dimP()-state_.dimTheta(),0,3,3) = skew(-it->pose.block<3,1>(0,3)); 
        G.block(G.rows()-state_.dimTheta()-3,0,3,3) = Eigen::Matrix3d::Identity();
      }
      P_aug = (F*P_aug*F.transpose() + G*it->covariance.block<3,3>(3,3)*G.transpose()).eval(); 

      // Update state and covariance
      state_.setX(X_aug); // TODO: move outside of loop (need to make loop independent of state_)
      state_.setP(P_aug);

      // Add to list of estimated contact positions
      estimated_contact_positions_.insert(pair<int,int> (it->id, startIndex));
    }
  }
}


// Create Observation from vector of landmark measurements
void InEKF::CorrectLandmarks(const vectorLandmarks& measured_landmarks) {
  Eigen::VectorXd Z, Y, b;
  Eigen::MatrixXd H, N, PI;
  vectorLandmarks new_landmarks;
  vector<int> used_landmark_ids;

  for (vectorLandmarksIterator it=measured_landmarks.begin(); it!=measured_landmarks.end(); ++it) {
    // Detect and skip if an ID is not unique (this would cause singularity issues in InEKF::Correct)
    if (find(used_landmark_ids.begin(), used_landmark_ids.end(), it->id) != used_landmark_ids.end()) { 
      cout << "Duplicate landmark ID detected! Skipping measurement.\n";
      continue; 
    } 
    else { 
      used_landmark_ids.push_back(it->id); 
    }
    // See if we can find id in prior_landmarks or estimated_landmarks
    mapIntVector3dIterator it_prior = prior_landmarks_.find(it->id);
    map<int,int>::iterator it_estimated = estimated_landmarks_.find(it->id);
    if (it_prior!=prior_landmarks_.end()) {
      // Found in prior landmark set
      const int dimX = state_.dimX();
      const int dimTheta = state_.dimTheta();
      const int dimP = state_.dimP();
      int startIndex;

      // Fill out H
      startIndex = H.rows();
      H.conservativeResize(startIndex+3, dimP);
      H.block(startIndex,0,3,dimP).setZero();
      if (state_.getStateType() == StateType::WorldCentric) {
          H.block(startIndex,0,3,3) = skew(it_prior->second); // skew(p_wl)
          H.block(startIndex,6,3,3) = -Eigen::Matrix3d::Identity(); // -I    
      } 
      else {
          H.block(startIndex,0,3,3) = skew(-it_prior->second); // -skew(p_wl)
          H.block(startIndex,6,3,3) = Eigen::Matrix3d::Identity(); // I    
      }

      // Fill out N
      startIndex = N.rows();
      N.conservativeResize(startIndex+3, startIndex+3);
      N.block(startIndex,0,3,startIndex).setZero();
      N.block(0,startIndex,startIndex,3).setZero();
      N.block(startIndex,startIndex,3,3) = state_.getWorldRotation() * it->covariance * state_.getWorldRotation().transpose();

      // Fill out Z
      startIndex = Z.rows();
      Z.conservativeResize(startIndex+3, Eigen::NoChange);
      const auto& R = state_.getRotation();
      const auto& p = state_.getPosition();
      const auto& l = state_.getVector(it_estimated->second);  
      if (state_.getStateType() == StateType::WorldCentric) {
          Z.segment(startIndex,3) = R*it->position - (l - it_prior->second); 
      } 
      else {
          Z.segment(startIndex,3) = R.transpose()*(it->position - (p - it_prior->second)); 
      }
    } 
    else if (it_estimated!=estimated_landmarks_.end()) {;
      // Found in estimated landmark set
      const int dimX = state_.dimX();
      const int dimTheta = state_.dimTheta();
      const int dimP = state_.dimP();
      int startIndex;

      // Fill out H
      startIndex = H.rows();
      H.conservativeResize(startIndex+3, dimP);
      H.block(startIndex,0,3,dimP).setZero();
      if (state_.getStateType() == StateType::WorldCentric) {
          H.block(startIndex,6,3,3) = -Eigen::Matrix3d::Identity(); // -I
          H.block(startIndex,3*it_estimated->second-dimTheta,3,3) = Eigen::Matrix3d::Identity(); // I
      } 
      else {
          H.block(startIndex,6,3,3) = Eigen::Matrix3d::Identity(); // I
          H.block(startIndex,3*it_estimated->second-dimTheta,3,3) = -Eigen::Matrix3d::Identity(); // -I
      }

      // Fill out N
      startIndex = N.rows();
      N.conservativeResize(startIndex+3, startIndex+3);
      N.block(startIndex,0,3,startIndex).setZero();
      N.block(0,startIndex,startIndex,3).setZero();
      N.block(startIndex,startIndex,3,3).noalias() = state_.getWorldRotation() * it->covariance * state_.getWorldRotation().transpose();

      // Fill out Z
      startIndex = Z.rows();
      Z.conservativeResize(startIndex+3, Eigen::NoChange);
      const auto& R = state_.getRotation();
      const auto& p = state_.getPosition();
      const auto& l = state_.getVector(it_estimated->second);  
      if (state_.getStateType() == StateType::WorldCentric) {
          Z.segment(startIndex,3).noalias() = R*it->position - (l - p); 
      } 
      else {
          Z.segment(startIndex,3).noalias() = R.transpose()*(it->position - (p - l)); 
      }
    } 
    else {
      // First time landmark as been detected (add to list for later state augmentation)
      new_landmarks.push_back(*it);
    }
  }

  // Correct state using stacked observation
  if (Z.rows()>0) {
    if (state_.getStateType() == StateType::WorldCentric) {
      this->CorrectRightInvariant(Z,H,N);
    } 
    else {
      this->CorrectLeftInvariant(Z,H,N);
    }
  }

    // Augment state with newly detected landmarks
  if (new_landmarks.size() > 0) {
    Eigen::MatrixXd X_aug = state_.getX(); 
    Eigen::MatrixXd P_aug = state_.getP();
    for (vectorLandmarksIterator it=new_landmarks.begin(); it!=new_landmarks.end(); ++it) {
      // Initialize new landmark mean
      const int startIndex = X_aug.rows();
      X_aug.conservativeResize(startIndex+1, startIndex+1);
      X_aug.block(startIndex,0,1,startIndex).setZero();
      X_aug.block(0,startIndex,startIndex,1).setZero();
      X_aug(startIndex, startIndex) = 1;
      X_aug.block(0,startIndex,3,1) = state_.getPosition() + state_.getRotation()*it->position;

      // Initialize new landmark covariance - TODO:speed up
      Eigen::MatrixXd F = Eigen::MatrixXd::Zero(state_.dimP()+3,state_.dimP()); 
      F.block(0,0,state_.dimP()-state_.dimTheta(),state_.dimP()-state_.dimTheta()).setIdentity(); // for old X
      F.block(state_.dimP()-state_.dimTheta()+3,state_.dimP()-state_.dimTheta(),state_.dimTheta(),state_.dimTheta()).setIdentity(); // for theta
      Eigen::MatrixXd G = Eigen::MatrixXd::Zero(F.rows(),3);
      // Blocks for new landmark
      if (error_type_==ErrorType::RightInvariant) {
        F.block(state_.dimP()-state_.dimTheta(),6,3,3) = Eigen::Matrix3d::Identity(); 
        G.block(G.rows()-state_.dimTheta()-3,0,3,3) = state_.getRotation();
      } else {
        F.block(state_.dimP()-state_.dimTheta(),6,3,3) = Eigen::Matrix3d::Identity(); 
        F.block(state_.dimP()-state_.dimTheta(),0,3,3) = skew(-it->position); 
        G.block(G.rows()-state_.dimTheta()-3,0,3,3) = Eigen::Matrix3d::Identity();
      }
      P_aug = (F*P_aug*F.transpose() + G*it->covariance*G.transpose()).eval();

      // Update state and covariance
      state_.setX(X_aug);
      state_.setP(P_aug);

      // Add to list of estimated landmarks
      estimated_landmarks_.insert(pair<int,int> (it->id, startIndex));
    }
  }
}


// Remove landmarks by IDs
void InEKF::RemoveLandmarks(const int landmark_id) {
    // Search for landmark in state
  map<int,int>::iterator it = estimated_landmarks_.find(landmark_id);
  if (it!=estimated_landmarks_.end()) {
    // Get current X and P
    Eigen::MatrixXd X_rem = state_.getX(); 
    Eigen::MatrixXd P_rem = state_.getP();
    // Remove row and column from X
    removeRowAndColumn(X_rem, it->second);
    // Remove 3 rows and columns from P
    int startIndex = 3 + 3*(it->second-3);
    removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
    removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
    removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
    // Update all indices for estimated_landmarks and estimated_contact_positions (TODO: speed this up)
    for (map<int,int>::iterator it2=estimated_landmarks_.begin(); it2!=estimated_landmarks_.end(); ++it2) {
      if (it2->second > it->second) it2->second -= 1;
    }
    for (map<int,int>::iterator it2=estimated_contact_positions_.begin(); it2!=estimated_contact_positions_.end(); ++it2) {
      if (it2->second > it->second) it2->second -= 1;
    }
    // Remove from list of estimated landmark positions (after we are done with iterator)
    estimated_landmarks_.erase(it->first);
    // Update state and covariance
    state_.setX(X_rem);
    state_.setP(P_rem);   
  }
}


// Remove landmarks by IDs
void InEKF::RemoveLandmarks(const std::vector<int>& landmark_ids) {
  // Loop over landmark_ids and remove
  for (int i=0; i<landmark_ids.size(); ++i) {
    this->RemoveLandmarks(landmark_ids[i]);
  }
}


// Keep landmarks by IDs
void InEKF::KeepLandmarks(const std::vector<int>& landmark_ids) {
  std::cout << std::endl;
  // Loop through estimated landmarks removing ones not found in the list
  std::vector<int> ids_to_erase;
  for (map<int,int>::iterator it=estimated_landmarks_.begin(); it!=estimated_landmarks_.end(); ++it) {
    std::vector<int>::const_iterator it_found = find(landmark_ids.begin(), landmark_ids.end(), it->first);
    if (it_found==landmark_ids.end()) {
      // Get current X and P
      Eigen::MatrixXd X_rem = state_.getX(); 
      Eigen::MatrixXd P_rem = state_.getP();
      // Remove row and column from X
      removeRowAndColumn(X_rem, it->second);
      // Remove 3 rows and columns from P
      int startIndex = 3 + 3*(it->second-3);
      removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
      removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
      removeRowAndColumn(P_rem, startIndex); // TODO: Make more efficient
      // Update all indices for estimated_landmarks and estimated_contact_positions (TODO: speed this up)
      for (map<int,int>::iterator it2=estimated_landmarks_.begin(); it2!=estimated_landmarks_.end(); ++it2) {
        if (it2->second > it->second) it2->second -= 1;
      }
      for (map<int,int>::iterator it2=estimated_contact_positions_.begin(); it2!=estimated_contact_positions_.end(); ++it2) {
        if (it2->second > it->second) it2->second -= 1;
      }
      // Add to list of ids to erase
      ids_to_erase.push_back(it->first);
      // Update state and covariance
      state_.setX(X_rem);
      state_.setP(P_rem);   
    }
  }
  // Remove from list of estimated landmark positions (after we are done with iterator)
  for (int i=0; i<ids_to_erase.size(); ++i) {
    estimated_landmarks_.erase(ids_to_erase[i]);
  }
}


// Remove prior landmarks by IDs
void InEKF::RemovePriorLandmarks(const int landmark_id) {
  // Search for landmark in state
  mapIntVector3dIterator it = prior_landmarks_.find(landmark_id);
  if (it!=prior_landmarks_.end()) { 
    // Remove from list of estimated landmark positions
    prior_landmarks_.erase(it->first);
  }
}


// Remove prior landmarks by IDs
void InEKF::RemovePriorLandmarks(const std::vector<int>& landmark_ids) {
  // Loop over landmark_ids and remove
  for (int i=0; i<landmark_ids.size(); ++i) {
    this->RemovePriorLandmarks(landmark_ids[i]);
  }
}


// Corrects state using magnetometer measurements (Right Invariant)
void InEKF::CorrectMagnetometer(const Eigen::Vector3d& measured_magnetic_field, const Eigen::Matrix3d& covariance) {
    // Eigen::VectorXd Y, b;
    // Eigen::MatrixXd H, N, PI;

    // // Get Rotation Estimate
    // Eigen::Matrix3d R = state_.getRotation();

    // // Fill out observation data
    // int dimX = state_.dimX();
    // int dimTheta = state_.dimTheta();
    // int dimP = state_.dimP();

    // // Fill out Y
    // Y.conservativeResize(dimX, Eigen::NoChange);
    // Y.segment(0,dimX) = Eigen::VectorXd::Zero(dimX);
    // Y.segment<3>(0) = measured_magnetic_field;

    // // Fill out b
    // b.conservativeResize(dimX, Eigen::NoChange);
    // b.segment(0,dimX) = Eigen::VectorXd::Zero(dimX);
    // b.segment<3>(0) = magnetic_field_;

    // // Fill out H
    // H.conservativeResize(3, dimP);
    // H.block(0,0,3,dimP) = Eigen::MatrixXd::Zero(3,dimP);
    // H.block<3,3>(0,0) = skew(magnetic_field_); 

    // // Fill out N
    // N.conservativeResize(3, 3);
    // N = R * covariance * R.transpose();

    // // Fill out PI      
    // PI.conservativeResize(3, dimX);
    // PI.block(0,0,3,dimX) = Eigen::MatrixXd::Zero(3,dimX);
    // PI.block(0,0,3,3) = Eigen::Matrix3d::Identity();
    

    // // Correct state using stacked observation
    // Observation obs(Y,b,H,N,PI);
    // if (!obs.empty()) {
    //     this->CorrectRightInvariant(obs);
    //     // cout << obs << endl;
    // }
}


// Observation of absolute position - GPS (Left-Invariant Measurement)
void InEKF::CorrectPosition(const Eigen::Vector3d& measured_position, 
                            const Eigen::Matrix3d& covariance, 
                            const Eigen::Vector3d& indices) {
  // Eigen::VectorXd Y, b;
  // Eigen::MatrixXd H, N, PI;

  // // Fill out observation data
  // int dimX = state_.dimX();
  // int dimTheta = state_.dimTheta();
  // int dimP = state_.dimP();

  // // Fill out Y
  // Y.conservativeResize(dimX, Eigen::NoChange);
  // Y.segment(0,dimX) = Eigen::VectorXd::Zero(dimX);
  // Y.segment<3>(0) = measured_position;
  // Y(4) = 1;       

  // // Fill out b
  // b.conservativeResize(dimX, Eigen::NoChange);
  // b.segment(0,dimX) = Eigen::VectorXd::Zero(dimX);
  // b(4) = 1;       

  // // Fill out H
  // H.conservativeResize(3, dimP);
  // H.block(0,0,3,dimP) = Eigen::MatrixXd::Zero(3,dimP);
  // H.block<3,3>(0,6) = Eigen::Matrix3d::Identity(); 

  // // Fill out N
  // N.conservativeResize(3, 3);
  // N = covariance;

  // // Fill out PI      
  // PI.conservativeResize(3, dimX);
  // PI.block(0,0,3,dimX) = Eigen::MatrixXd::Zero(3,dimX);
  // PI.block(0,0,3,3) = Eigen::Matrix3d::Identity();

  // // Modify measurement based on chosen indices
  // const double HIGH_UNCERTAINTY = 1e6;
  // Eigen::Vector3d p = state_.getPosition();
  // if (!indices(0)) { 
  //   Y(0) = p(0);
  //   N(0,0) = HIGH_UNCERTAINTY;
  //   N(0,1) = 0;
  //   N(0,2) = 0;
  //   N(1,0) = 0;
  //   N(2,0) = 0;
  //   } 
  // if (!indices(1)) { 
  //   Y(1) = p(1);
  //   N(1,0) = 0;
  //   N(1,1) = HIGH_UNCERTAINTY;
  //   N(1,2) = 0;
  //   N(0,1) = 0;
  //   N(2,1) = 0;
  //   } 
  // if (!indices(2)) { 
  //   Y(2) = p(2);
  //   N(2,0) = 0;
  //   N(2,1) = 0;
  //   N(2,2) = HIGH_UNCERTAINTY;
  //   N(0,2) = 0;
  //   N(1,2) = 0;
  //   } 

  // // Correct state using stacked observation
  // Observation obs(Y,b,H,N,PI);
  // if (!obs.empty()) {
  //   this->CorrectLeftInvariant(obs);
  //   // cout << obs << endl;
  // }
}


// Observation of absolute z-position of contact points (Left-Invariant Measurement)
void InEKF::CorrectContactPosition(const int id, 
                                   const Eigen::Vector3d& measured_contact_position, 
                                   const Eigen::Matrix3d& covariance, 
                                   const Eigen::Vector3d& indices) {
  Eigen::VectorXd Z_full, Z;
  Eigen::MatrixXd PI, H_full, N_full, H, N;

  // See if we can find id estimated_contact_positions
  map<int,int>::iterator it_estimated = estimated_contact_positions_.find(id);
  if (it_estimated!=estimated_contact_positions_.end()) { 

    // Fill out PI
    int startIndex;
    for (int i=0; i<3; ++i) {
      if (indices(i) != 0) {
        startIndex = PI.rows();
        PI.conservativeResize(startIndex+1, 3);  
        PI.template block<1,3>(startIndex,0).setZero();
        PI.coeffRef(startIndex,i) = 1;
      }  
    }
    if (PI.rows()==0) {
      return;
    }

    // Fill out observation data
    const int dimX = state_.dimX();
    const int dimTheta = state_.dimTheta();
    const int dimP = state_.dimP();

    // Get contact position
    const auto& d = state_.getVector(it_estimated->second);

    // Fill out H
    H_full = Eigen::MatrixXd::Zero(3,dimP);
    H_full.block<3,3>(0,0) = -skew(d);
    H_full.block<3,3>(0,3*it_estimated->second-6) = Eigen::Matrix3d::Identity();
    H.noalias() = PI*H_full;

    // Fill out N
    N_full = covariance;   
    N.noalias() = PI*N_full*PI.transpose();

    // Fill out Z
    Z_full = measured_contact_position - d; 
    Z.noalias() = PI*Z_full;

    // Correct
    this->CorrectRightInvariant(Z,H,N);
  }
}


void removeRowAndColumn(Eigen::MatrixXd& M, int index) {
  const unsigned int dimX = M.cols();
  // cout << "Removing index: " << index<< endl;
  M.block(index,0,dimX-index-1,dimX) = M.bottomRows(dimX-index-1).eval();
  M.block(0,index,dimX,dimX-index-1) = M.rightCols(dimX-index-1).eval();
  M.conservativeResize(dimX-1,dimX-1);
}

} // end inekf namespace
