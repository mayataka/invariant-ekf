#ifndef INEKF_DISCRETE_NOISE_MATRIX_HPP_
#define INEKF_DISCRETE_NOISE_MATRIX_HPP_

#include <map>

#include "Eigen/Core"

#include "inekf/inekf_state.hpp"
#include "inekf/noise_params.hpp"
#include "inekf/lie_group.hpp"
#include "inekf/macros.hpp"
#include "inekf/error_type.hpp"


namespace inekf {

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

} // namespace inekf 

#endif // INEKF_DISCRETE_NOISE_MATRIX_HPP_