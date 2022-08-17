#include "legged_state_estimator/inekf/state_transition_matrix.hpp"


namespace legged_state_estimator {

StateTransitionMatrix::StateTransitionMatrix() { 
  phi_.setZero();
  G0_.setZero();
  G1_.setZero();
  G2_.setZero();
  G0t_.setZero();
  G1t_.setZero();
  G2t_.setZero();
  G3t_.setZero();
  ax_.setZero();
  wx_.setZero();
  wx2_.setZero();
  Phi25L_.setZero();
  Phi35L_.setZero();
  g_ << 0, 0, -9.81;
  gx_.setZero();
  RG0_.setZero();
  RG1dt_.setZero();
  RG2dt2_.setZero();
}


StateTransitionMatrix::StateTransitionMatrix(const InEKFState& state)
  : StateTransitionMatrix() {
  const int dimP = state.dimP();
  Phi_.resize(dimP, dimP);
  Phi_.setZero();
}


StateTransitionMatrix::StateTransitionMatrix(const InEKFState& state, 
                                             const ErrorType error_type)
  : StateTransitionMatrix(state) {
  error_type_ = error_type;
}


void StateTransitionMatrix::compute(const InEKFState& state, const Eigen::Vector3d& w, 
                                    const Eigen::Vector3d& a, double dt) {
  phi_ = w*dt;
  Gamma_SO3(phi_, G0_, 0); // Computation can be sped up by computing G0,G1,G2 all at once
  Gamma_SO3(phi_, G1_, 1); // TODO: These are also needed for the mean propagation, we should not compute twice
  Gamma_SO3(phi_, G2_, 2);
  G0t_ = G0_.transpose();
  G1t_ = G1_.transpose();
  G2t_ = G2_.transpose();
  Gamma_SO3(-phi_, G3t_, 3);

  // Compute the complicated bias terms (derived for the left invariant case)
  skew(a, ax_);
  skew(w, wx_);
  wx2_.noalias() = wx_ * wx_;
  const double dt2 = dt*dt;
  const double dt3 = dt2*dt;
  const double theta = w.norm();
  const double theta2 = theta*theta;
  const double theta3 = theta2*theta;
  const double theta4 = theta3*theta;
  const double theta5 = theta4*theta;
  const double theta6 = theta5*theta;
  const double theta7 = theta6*theta;
  const double thetadt = theta*dt;
  const double thetadt2 = thetadt*thetadt;
  const double thetadt3 = thetadt2*thetadt;
  const double sinthetadt = std::sin(thetadt);
  const double costhetadt = std::cos(thetadt);
  const double sin2thetadt = std::sin(2*thetadt);
  const double cos2thetadt = std::cos(2*thetadt);
  const double thetadtcosthetadt = thetadt*costhetadt;
  const double thetadtsinthetadt = thetadt*sinthetadt;

  Phi25L_.noalias() = G0t_*(ax_*G2t_*dt2 
      + ((sinthetadt-thetadtcosthetadt)/(theta3))*(wx_*ax_)
      - ((cos2thetadt-4*costhetadt+3)/(4*theta4))*(wx_*ax_*wx_)
      + ((4*sinthetadt+sin2thetadt-4*thetadtcosthetadt-2*thetadt)/(4*theta5))*(wx_*ax_*wx2_)
      + ((thetadt2-2*thetadtsinthetadt-2*costhetadt+2)/(2*theta4))*(wx2_*ax_)
      - ((6*thetadt-8*sinthetadt+sin2thetadt)/(4*theta5))*(wx2_*ax_*wx_)
      + ((2*thetadt2-4*thetadtsinthetadt-cos2thetadt+1)/(4*theta6))*(wx2_*ax_*wx2_) );
  Phi35L_.noalias() = G0t_*(ax_*G3t_*dt3
      - ((thetadtsinthetadt+2*costhetadt-2)/(theta4))*(wx_*ax_)
      - ((6*thetadt-8*sinthetadt+sin2thetadt)/(8*theta5))*(wx_*ax_*wx_)
      - ((2*thetadt2+8*thetadtsinthetadt+16*costhetadt+cos2thetadt-17)/(8*theta6))*(wx_*ax_*wx2_)
      + ((thetadt3+6*thetadt-12*sinthetadt+6*thetadtcosthetadt)/(6*theta5))*(wx2_*ax_)
      - ((6*thetadt2+16*costhetadt-cos2thetadt-15)/(8*theta6))*(wx2_*ax_*wx_)
      + ((4*thetadt3+6*thetadt-24*sinthetadt-3*sin2thetadt+24*thetadtcosthetadt)/(24*theta7))*(wx2_*ax_*wx2_) );

  // TODO: Get better approximation using taylor series when theta < tol
  const double tol =  1e-6;
  if (theta < tol) {
    Phi25L_ = (1.0/2.0)*ax_*dt2;
    Phi35L_ = (1.0/6.0)*ax_*dt3;
  }

  // Fill out analytical state transition matrices
  const int dimX = state.dimX();
  const int dimTheta = state.dimTheta();
  const int dimP = state.dimP();
  if (Phi_.cols() != dimP || Phi_.rows() != dimP) {
    Phi_.resize(dimP, dimP);
  }
  Phi_.setZero();
  if  ((state.getStateType() == StateType::WorldCentric && error_type_ == ErrorType::LeftInvariant) || 
        (state.getStateType() == StateType::BodyCentric && error_type_ == ErrorType::RightInvariant)) {
    // Compute left-invariant state transisition matrix
    Phi_.template block<3,3>(0,0) = G0t_; // Phi_11
    Phi_.template block<3,3>(3,0).noalias() = -G0t_ * skew(G1_*a) * dt; // Phi_21
    Phi_.template block<3,3>(6,0).noalias() = -G0t_ * skew(G2_*a) * dt2; // Phi_31
    Phi_.template block<3,3>(3,3) = G0t_; // Phi_22
    Phi_.template block<3,3>(6,3) = G0t_*dt; // Phi_32
    Phi_.template block<3,3>(6,6) = G0t_; // Phi_33
    for (int i=5; i<dimX; ++i) {
      Phi_.template block<3,3>((i-2)*3,(i-2)*3) = G0t_; // Phi_(3+i)(3+i)
    }
    Phi_.template block<3,3>(0,dimP-dimTheta) = -G1t_ * dt; // Phi_15
    Phi_.template block<3,3>(3,dimP-dimTheta) = Phi25L_; // Phi_25
    Phi_.template block<3,3>(6,dimP-dimTheta) = Phi35L_; // Phi_35
    Phi_.template block<3,3>(3,dimP-dimTheta+3) = -G1t_ * dt; // Phi_26
    Phi_.template block<3,3>(6,dimP-dimTheta+3).noalias() = -G0t_ * G2_ * dt2; // Phi_36
  } 
  else {
    // Compute right-invariant state transition matrix (Assumes unpropagated state)
    skew(g_, gx_);
    const auto& R = state.getRotation();
    const auto& v = state.getVelocity();
    const auto& p = state.getPosition();
    RG0_.noalias() = R*G0_;
    RG1dt_.noalias() = R*G1_*dt;
    RG2dt2_.noalias() = R*G2_*dt2;
    Phi_.template block<3,3>(3,0) = gx_*dt; // Phi_21
    Phi_.template block<3,3>(6,0) = 0.5*gx_*dt2; // Phi_31
    Phi_.template block<3,3>(6,3) = Eigen::Matrix3d::Identity()*dt; // Phi_32
    Phi_.template block<3,3>(0,dimP-dimTheta) = -RG1dt_; // Phi_15
    Phi_.template block<3,3>(3,dimP-dimTheta).noalias() = -skew(v+RG1dt_*a+g_*dt)*RG1dt_ + RG0_*Phi25L_; // Phi_25
    Phi_.template block<3,3>(6,dimP-dimTheta).noalias() = -skew(p+v*dt+RG2dt2_*a+0.5*g_*dt2)*RG1dt_ + RG0_*Phi35L_; // Phi_35
    for (int i=5; i<dimX; ++i) {
      Phi_.template block<3,3>((i-2)*3,dimP-dimTheta).noalias() = -skew(state.getVector(i))*RG1dt_; // Phi_(3+i)5
    }
    Phi_.template block<3,3>(3,dimP-dimTheta+3) = -RG1dt_; // Phi_26
    Phi_.template block<3,3>(6,dimP-dimTheta+3) = -RG2dt2_; // Phi_36
  }
}

} // namespace legged_state_estimator 
