#ifndef INEKF_DISCRETE_NOISE_MATRIX_HPP_
#define INEKF_DISCRETE_NOISE_MATRIX_HPP_

#include <map>

#include "Eigen/Core"

#include "legged_state_estimator/inekf/inekf_state.hpp"
#include "legged_state_estimator/inekf/noise_params.hpp"
#include "legged_state_estimator/inekf/lie_group.hpp"
#include "legged_state_estimator/inekf/error_type.hpp"
#include "legged_state_estimator/macros.hpp"


namespace legged_state_estimator {

class DiscreteNoiseMatrix {
public:
  DiscreteNoiseMatrix();

  DiscreteNoiseMatrix(const InEKFState& state);

  DiscreteNoiseMatrix(const InEKFState& state, const ErrorType error_type);

  INEKF_USE_DEFAULT_DESTTUCTOR(DiscreteNoiseMatrix);
  INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(DiscreteNoiseMatrix);
  INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(DiscreteNoiseMatrix);
  INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(DiscreteNoiseMatrix);
  INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(DiscreteNoiseMatrix);

  void compute(const InEKFState& state, const NoiseParams& noise_params, 
               const std::map<int,int>& estimated_contact_positions, 
               const Eigen::MatrixXd& Phi, double dt);

  const Eigen::MatrixXd& Qd() const {
    return Qd_;
  }

private:
  ErrorType error_type_ = ErrorType::LeftInvariant; 
  Eigen::MatrixXd G_, Qc_, PhiG_, Qd_;

};

} // namespace legged_state_estimator 

#endif // INEKF_DISCRETE_NOISE_MATRIX_HPP_