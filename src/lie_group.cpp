/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley <m.ross.hartley@gmail.com>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   LieGroup.cpp
 *  @author Ross Hartley
 *  @brief  Source file for various Lie Group functions 
 *  @date   September 25, 2018
 **/

#include "legged_state_estimator/inekf/lie_group.hpp"


namespace legged_state_estimator {

using namespace std;

long int factorial(const int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
  // Convert vector to skew-symmetric matrix
  Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
  skew(v, M);
  return M;
}

void skew(const Eigen::Vector3d& v, Eigen::Matrix3d& M) {
  // Convert vector to skew-symmetric matrix
  M <<           0, -v.coeff(2),  v.coeff(1),
        v.coeff(2),           0, -v.coeff(0), 
       -v.coeff(1),  v.coeff(0),           0;
}


Eigen::Matrix3d Gamma_SO3(const Eigen::Vector3d& w, const int m, 
                          const double exp_map_tol) {
  Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
  Gamma_SO3(w, M, m, exp_map_tol);
  return M;
}


void Gamma_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M, const int m, 
               const double exp_map_tol) {
  // Computes mth integral of the exponential map: \Gamma_m = \sum_{n=0}^{\infty} \dfrac{1}{(n+m)!} (w^\wedge)^n
  assert(m>=0);
  const double theta = w.norm();
  if (theta < exp_map_tol) {
    M = (1.0/factorial(m))*Eigen::Matrix3d::Identity(); // TODO: There is a better small value approximation for exp() given in Trawny p.19
    return;
  } 
  const Eigen::Matrix3d A = skew(w);
  const double theta2 =  theta*theta;

  // Closed form solution for the first 3 cases
  switch (m) {
    case 0: // Exp map of SO(3)
      M = Eigen::Matrix3d::Identity() + (sin(theta)/theta)*A + ((1-cos(theta))/theta2)*A*A;
      return;
    
    case 1: // Left Jacobian of SO(3)
      // eye(3) - A*(1/theta^2) * (R - eye(3) - A);
      // eye(3) + (1-cos(theta))/theta^2 * A + (theta-sin(theta))/theta^3 * A^2;
      M = Eigen::Matrix3d::Identity() + ((1-cos(theta))/theta2)*A + ((theta-sin(theta))/(theta2*theta))*A*A;
      return;

    case 2: 
      // 0.5*eye(3) - (1/theta^2) * (R - eye(3) - A - 0.5*A^2);
      // 0.5*eye(3) + (theta-sin(theta))/theta^3 * A + (2*(cos(theta)-1) + theta^2)/(2*theta^4) * A^2
      M = 0.5*Eigen::Matrix3d::Identity() + (theta-sin(theta))/(theta2*theta)*A + (theta2 + 2*cos(theta)-2)/(2*theta2*theta2)*A*A;
      return;

    default: // General case 
      const Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + (sin(theta)/theta)*A + ((1-cos(theta))/theta2)*A*A;
      Eigen::Matrix3d S = Eigen::Matrix3d::Identity();
      Eigen::Matrix3d Ak = Eigen::Matrix3d::Identity();
      long int kfactorial = 1;
      for (int k=1; k<=m; ++k) {
        kfactorial = kfactorial*k;
        Ak = (Ak*A).eval();
        S = (S + (1.0/kfactorial)*Ak).eval();
      }
      if (m==0) { 
          M = R;
          return;
      } 
      else if (m%2){ // odd 
          M = (1.0/kfactorial)*Eigen::Matrix3d::Identity() + (pow(-1,(m+1)/2)/pow(theta,m+1))*A * (R - S);
          return;
      } 
      else { // even
          M = (1.0/kfactorial)*Eigen::Matrix3d::Identity() + (pow(-1,m/2)/pow(theta,m)) * (R - S);
          return;
      }
  }
}

Eigen::Matrix3d Exp_SO3(const Eigen::Vector3d& w) {
  // Computes the vectorized exponential map for SO(3)
  return Gamma_SO3(w, 0);
}

void Exp_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M) {
  // Computes the vectorized exponential map for SO(3)
  Gamma_SO3(w, M, 0);
}

Eigen::Matrix3d LeftJacobian_SO3(const Eigen::Vector3d& w) {
  // Computes the Left Jacobian of SO(3)
  return Gamma_SO3(w, 1);
}

void LeftJacobian_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M) {
  // Computes the Left Jacobian of SO(3)
  Gamma_SO3(w, M, 1);
}

Eigen::Matrix3d RightJacobian_SO3(const Eigen::Vector3d& w) {
  // Computes the Right Jacobian of SO(3)
  return Gamma_SO3(-w, 1);
}

void RightJacobian_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M) {
  // Computes the Right Jacobian of SO(3)
  Gamma_SO3(-w, M, 1);
}

Eigen::MatrixXd Exp_SEK3(const Eigen::VectorXd& v,  const double exp_map_tol) {
  // Computes the vectorized exponential map for SE_K(3)
  const int K = (v.size()-3)/3;
  Eigen::MatrixXd X = Eigen::MatrixXd::Identity(3+K,3+K);
  Eigen::Matrix3d R;
  Eigen::Matrix3d Jl;
  const auto& w = v.head(3);
  double theta = w.norm();
  if (theta < exp_map_tol) {
    R = Eigen::Matrix3d::Identity();
    Jl = Eigen::Matrix3d::Identity();
  } else {
    const Eigen::Matrix3d A = skew(w);
    const double theta2 = theta*theta;
    const double stheta = sin(theta);
    const double ctheta = cos(theta);
    const double oneMinusCosTheta2 = (1-ctheta)/(theta2);
    const Eigen::Matrix3d A2 = A*A;
    R.noalias() = Eigen::Matrix3d::Identity() 
                   + (stheta/theta) * A + oneMinusCosTheta2 * A2;
    Jl.noalias() = Eigen::Matrix3d::Identity() 
                   + oneMinusCosTheta2*A + ((theta-stheta)/(theta2*theta)) * A2;
  }
  X.block<3,3>(0,0) = R;
  for (int i=0; i<K; ++i) {
      X.block<3,1>(0,3+i).noalias() = Jl * v.segment<3>(3+3*i);
  }
  return X;
}

void Exp_SEK3(const Eigen::VectorXd& v, Eigen::MatrixXd& X, const double exp_map_tol) {
  // Computes the vectorized exponential map for SE_K(3)
  const int K = (v.size()-3)/3;
  if (X.cols() != K || X.rows() != K) {
    X.resize(3+K, 3+K);
  }
  X.setIdentity();
  Eigen::Matrix3d R;
  Eigen::Matrix3d Jl;
  const auto& w = v.head(3);
  double theta = w.norm();
  if (theta < exp_map_tol) {
    R = Eigen::Matrix3d::Identity();
    Jl = Eigen::Matrix3d::Identity();
  } else {
    const Eigen::Matrix3d A = skew(w);
    const double theta2 = theta*theta;
    const double stheta = sin(theta);
    const double ctheta = cos(theta);
    const double oneMinusCosTheta2 = (1-ctheta)/(theta2);
    const Eigen::Matrix3d A2 = A*A;
    R.noalias() = Eigen::Matrix3d::Identity() 
                   + (stheta/theta) * A + oneMinusCosTheta2 * A2;
    Jl.noalias() = Eigen::Matrix3d::Identity() 
                   + oneMinusCosTheta2*A + ((theta-stheta)/(theta2*theta)) * A2;
  }
  X.block<3,3>(0,0) = R;
  for (int i=0; i<K; ++i) {
      X.block<3,1>(0,3+i).noalias() = Jl * v.segment<3>(3+3*i);
  }
}

Eigen::MatrixXd Adjoint_SEK3(const Eigen::MatrixXd& X) {
  // Compute Adjoint(X) for X in SE_K(3)
  const int K = X.cols()-3;
  Eigen::MatrixXd Adj = Eigen::MatrixXd::Zero(3+3*K, 3+3*K);
  const auto& R = X.block<3,3>(0,0);
  Adj.block<3,3>(0,0) = R;
  for (int i=0; i<K; ++i) {
    Adj.block<3,3>(3+3*i,3+3*i) = R;
    Adj.block<3,3>(3+3*i,0).noalias() = skew(X.block<3,1>(0,3+i)) * R;
  }
  return Adj;
}

void Adjoint_SEK3(const Eigen::MatrixXd& X, Eigen::MatrixXd& Adj) {
  // Compute Adjoint(X) for X in SE_K(3)
  const int K = X.cols()-3;
  if (Adj.cols() != 3+3*K || Adj.rows() != 3+3*K) {
    Adj.setZero(3+3*K, 3+3*K);
  }
  const auto& R = X.block<3,3>(0,0);
  Adj.block<3,3>(0,0) = R;
  for (int i=0; i<K; ++i) {
    Adj.block<3,3>(3+3*i,3+3*i) = R;
    Adj.block<3,3>(3+3*i,0).noalias() = skew(X.block<3,1>(0,3+i)) * R;
  }
}

} // namespace legged_state_estimator 