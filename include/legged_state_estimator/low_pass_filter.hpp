#ifndef LEGGED_STATE_ESTIMATOR_LOW_PASS_FILTER_HPP_
#define LEGGED_STATE_ESTIMATOR_LOW_PASS_FILTER_HPP_

#include <cmath>
#include <stdexcept>
#include <iostream>

#include "Eigen/Core"


namespace legged_state_estimator {

template <typename Scalar, int dim=Eigen::Dynamic>
class LowPassFilter {
public:
  using Vector = Eigen::Matrix<Scalar, dim, 1>;

  LowPassFilter(const Scalar sampling_time, const Scalar cutoff_freq,
                const int dynamic_size=0)
    : estimate_(),
      alpha_(0.0) {
    try {
      if (sampling_time <= 0) {
        throw std::out_of_range(
            "Invalid argment: sampling_time must be positive!");
      }
      if (cutoff_freq <= 0) {
        throw std::out_of_range(
            "Invalid argment: cutoff_freq must be positive!");
      }
      if (dim == Eigen::Dynamic && dynamic_size <= 0) {
        throw std::out_of_range(
            "Invalid argment: dynamic_size must be positive!");
      }
    }
    catch(const std::exception& e) {
      std::cerr << e.what() << '\n';
      std::exit(EXIT_FAILURE);
    }
    const Scalar tau = 1.0 / (2.0*M_PI*cutoff_freq);
    alpha_ = tau / (tau + sampling_time);
    if (dim == Eigen::Dynamic) {
      estimate_.resize(dynamic_size);
    }
    estimate_.setZero();
  }

  LowPassFilter()
    : estimate_(),
      alpha_(0.0) {
  }

  ~LowPassFilter() = default;

  void reset() {
    estimate_.setZero();
  }

  void reset(const Vector& estimate) {
    estimate_ = estimate;
  }

  void update(const Vector& obs) {
    estimate_.array() *= alpha_;
    estimate_.noalias() += (1.0-alpha_) * obs;
  }

  const Vector& getEstimate() const {
    return estimate_;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Vector estimate_;
  Scalar alpha_;
};

} // namespace legged_state_estimator 

#endif // LEGGED_STATE_ESTIMATOR_LOW_PASS_FILTER_HPP_ 