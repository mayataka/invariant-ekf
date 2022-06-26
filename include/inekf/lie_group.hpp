/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   lie_group.h
 *  @author Ross Hartley
 *  @brief  Header file for various Lie Group functions 
 *  @date   September 25, 2018
 **/

#ifndef INEKF_LIEGROUP_HPP_
#define INEKF_LIEGROUP_HPP_
#include <iostream>
#include <cmath>

#include "Eigen/Core"

namespace inekf {

long int factorial(const int n);
Eigen::Matrix3d skew(const Eigen::Vector3d& v);
void skew(const Eigen::Vector3d& v, Eigen::Matrix3d& M);

Eigen::Matrix3d Gamma_SO3(const Eigen::Vector3d& w, const int n,
                          const double exp_map_tol=1.0e-10);
void Gamma_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M, const int n,
               const double exp_map_tol=1.0e-10);

Eigen::Matrix3d Exp_SO3(const Eigen::Vector3d& w);
void Exp_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M);

Eigen::Matrix3d LeftJacobian_SO3(const Eigen::Vector3d& w);
void LeftJacobian_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M);

Eigen::Matrix3d RightJacobian_SO3(const Eigen::Vector3d& w);
void RightJacobian_SO3(const Eigen::Vector3d& w, Eigen::Matrix3d& M);

Eigen::MatrixXd Exp_SEK3(const Eigen::VectorXd& v, 
                         const double exp_map_tol=1.0e-10);
void Exp_SEK3(const Eigen::VectorXd& v, Eigen::Matrix3d& M, 
              const double exp_map_tol=1.0e-10);

Eigen::MatrixXd Adjoint_SEK3(const Eigen::MatrixXd& X);
void Adjoint_SEK3(const Eigen::MatrixXd& X, Eigen::MatrixXd& Adj);

} // namespace inekf 

#endif // INEKF_LIEGROUP_HPP_