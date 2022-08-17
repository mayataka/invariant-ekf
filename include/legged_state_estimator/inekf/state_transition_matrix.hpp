#ifndef INEKF_STATE_TRANSITION_MATRIX_HPP_
#define INEKF_STATE_TRANSITION_MATRIX_HPP_

#include "Eigen/Core"

#include "legged_state_estimator/inekf/inekf_state.hpp"
#include "legged_state_estimator/inekf/lie_group.hpp"
#include "legged_state_estimator/inekf/error_type.hpp"
#include "legged_state_estimator/macros.hpp"


namespace legged_state_estimator {

class StateTransitionMatrix {
public:
  StateTransitionMatrix();

  StateTransitionMatrix(const InEKFState& state);

  StateTransitionMatrix(const InEKFState& state, const ErrorType error_type);

  INEKF_USE_DEFAULT_DESTTUCTOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(StateTransitionMatrix);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(StateTransitionMatrix);

  void compute(const InEKFState& state, const Eigen::Vector3d& w, 
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

} // namespace legged_state_estimator 

#endif // INEKF_STATE_TRANSITION_MATRIX_HPP_