#ifndef INEKF_STATE_TRANSITION_MATRIX_HPP_
#define INEKF_STATE_TRANSITION_MATRIX_HPP_

#include "Eigen/Core"

#include "inekf/robot_state.hpp"
#include "inekf/lie_group.hpp"
#include "inekf/macros.hpp"
#include "inekf/error_type.hpp"


namespace inekf {

class StateTransitionMatrix {
public:
  StateTransitionMatrix();

  StateTransitionMatrix(const RobotState& state);

  StateTransitionMatrix(const RobotState& state, const ErrorType error_type);

  INEKF_USE_DEFAULT_DESTTUCTOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(StateTransitionMatrix);

  void compute(const RobotState& state, const Eigen::Vector3d& w, 
               const Eigen::Vector3d& a, double dt);

  const Eigen::MatrixXd& Phi() const {
    return Phi_;
  }

private:
  ErrorType error_type_ = ErrorType::LeftInvariant; 
  Eigen::Vector3d phi_;
  Eigen::Matrix3d G0_, G1_, G2_, G0t_, G1t_, G2t_, G3t_;
  Eigen::Matrix3d ax_, wx_, wx2_;
  Eigen::Matrix3d Phi25L_, Phi35L_;
  Eigen::Vector3d g_;
  Eigen::Matrix3d gx_;
  Eigen::Matrix3d RG0_, RG1dt_, RG2dt2_;
  Eigen::MatrixXd Phi_;

};

} // namespace inekf 

#endif // INEKF_STATE_TRANSITION_MATRIX_HPP_