#include "legged_state_estimator/inekf/discrete_noise_matrix.hpp"


namespace legged_state_estimator {

DiscreteNoiseMatrix::DiscreteNoiseMatrix() { 
}


DiscreteNoiseMatrix::DiscreteNoiseMatrix(const InEKFState& state)
  : DiscreteNoiseMatrix() {
  const int dimP = state.dimP();
  G_.resize(dimP, dimP);
  G_.setZero();
  Qc_.resize(dimP, dimP);
  Qc_.setZero();
  PhiG_.resize(dimP, dimP);
  PhiG_.setZero();
  Qd_.resize(dimP, dimP);
  Qd_.setZero();
}


DiscreteNoiseMatrix::DiscreteNoiseMatrix(const InEKFState& state, 
                                         const ErrorType error_type)
  : DiscreteNoiseMatrix(state) {
  error_type_ = error_type;
}


void DiscreteNoiseMatrix::compute(const InEKFState& state, const NoiseParams& noise_params,
                                  const std::map<int,int>& estimated_contact_positions, 
                                  const Eigen::MatrixXd& Phi, double dt) {
  const int dimX = state.dimX();
  const int dimTheta = state.dimTheta();
  const int dimP = state.dimP();    
  if (G_.rows() != dimP || G_.cols() != dimP) {
    G_.resize(dimP, dimP);
  }
  if (Qc_.rows() != dimP || Qc_.cols() != dimP) {
    Qc_.resize(dimP, dimP);
  }
  if (PhiG_.rows() != dimP || PhiG_.cols() != dimP) {
    PhiG_.resize(dimP, dimP);
  }
  if (Qd_.rows() != dimP || Qd_.cols() != dimP) {
    Qd_.resize(dimP, dimP);
  }

  G_.setIdentity();
  // Compute G using Adjoint of Xk if needed, otherwise identity (Assumes unpropagated state)
  if ((state.getStateType() == StateType::WorldCentric && error_type_ == ErrorType::RightInvariant) || 
      (state.getStateType() == StateType::BodyCentric && error_type_ == ErrorType::LeftInvariant)) {
    G_.block(0,0,dimP-dimTheta,dimP-dimTheta) = Adjoint_SEK3(state.getWorldX()); 
  }

  // Continuous noise covariance 
  Qc_.setZero();
  Qc_.template block<3,3>(0,0) = noise_params.getGyroscopeCov(); 
  Qc_.template block<3,3>(3,3) = noise_params.getAccelerometerCov();
  for (auto& e : estimated_contact_positions) {
    Qc_.template block<3,3>(3+3*(e.second-3),3+3*(e.second-3)) = noise_params.getContactCov(); // Contact noise terms
  }
  // TODO: Use kinematic orientation to map noise from contact frame to body frame (not needed if noise is isotropic)
  Qc_.template block<3,3>(dimP-dimTheta,dimP-dimTheta) = noise_params.getGyroscopeBiasCov();
  Qc_.template block<3,3>(dimP-dimTheta+3,dimP-dimTheta+3) = noise_params.getAccelerometerBiasCov();

  // Noise Covariance Discretization
  PhiG_.noalias() = Phi * G_;
  Qd_.noalias() = PhiG_ * Qc_ * PhiG_.transpose() * dt; // Approximated discretized noise matrix (TODO: compute analytical)
}

} // namespace legged_state_estimator 
