#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "inekf/slip_estimator.hpp"


namespace inekf {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(slip_estimator, m) {
  py::class_<SlipEstimatorSettings>(m, "SlipEstimatorSettings")
    .def(py::init<>())
    .def_readwrite("beta0", &SlipEstimatorSettings::beta0)
    .def_readwrite("beta1", &SlipEstimatorSettings::beta1)
    .def_readwrite("slip_velocity_cov_alpha", &SlipEstimatorSettings::slip_velocity_cov_alpha)
    .def_readwrite("slip_prob_threshold", &SlipEstimatorSettings::slip_prob_threshold)
    .def_readwrite("lpf_contact_surface_normal_cutoff", &SlipEstimatorSettings::lpf_contact_surface_normal_cutoff)
    .def_readwrite("lpf_friction_coefficient_cutoff", &SlipEstimatorSettings::lpf_friction_coefficient_cutoff);

  py::class_<SlipEstimator>(m, "SlipEstimator")
    .def(py::init<const RobotModel&, const SlipEstimatorSettings&, const double>(),
          py::arg("robot_model"), py::arg("settings"), py::arg("dt"))
    .def(py::init<>())
    .def("reset", &SlipEstimator::reset)
    .def("update", &SlipEstimator::update,
          py::arg("robot_model"), py::arg("contact_estimator"))
    .def("get_slip_state", &SlipEstimator::getSlipState)
    .def("get_slip_probability", &SlipEstimator::getSlipProbability)
    .def("get_slip_velocity_norm", &SlipEstimator::getSlipVelocityNorm)
    .def("get_slip_velocity_covariance", &SlipEstimator::getSlipVelocityCovariance)
    .def("get_friction_coefficient_estimate", &SlipEstimator::getFrictionCoefficientEstimate)
    .def("get_contact_surface_normal_estimate", &SlipEstimator::getContactSurfaceNormalEstimate)
    .def("get_contact_surface_estimate", &SlipEstimator::getContactSurfaceEstimate)
     .def("__str__", [](const SlipEstimator& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace inekf