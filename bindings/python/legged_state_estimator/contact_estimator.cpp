#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "legged_state_estimator/contact_estimator.hpp"


namespace legged_state_estimator {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_estimator, m) {
  py::class_<ContactEstimatorSettings>(m, "ContactEstimatorSettings")
    .def(py::init<>())
    .def_readwrite("beta0", &ContactEstimatorSettings::beta0)
    .def_readwrite("beta1", &ContactEstimatorSettings::beta1)
    .def_readwrite("contact_force_cov_alpha", &ContactEstimatorSettings::contact_force_cov_alpha)
    .def_readwrite("contact_prob_threshold", &ContactEstimatorSettings::contact_prob_threshold);

  py::class_<ContactEstimator>(m, "ContactEstimator")
    .def(py::init<const RobotModel&, const ContactEstimatorSettings&>(),
          py::arg("robot_model"), py::arg("settings"))
    .def(py::init<>())
    .def("reset", &ContactEstimator::reset)
    .def("update", &ContactEstimator::update,
          py::arg("robot_model"), py::arg("tauJ"))
    .def("get_contact_state", &ContactEstimator::getContactState)
    .def("get_contact_force_estimate", &ContactEstimator::getContactForceEstimate)
    .def("get_normal_contact_force_estimate", &ContactEstimator::getNormalContactForceEstimate)
    .def("get_contact_probability", &ContactEstimator::getContactProbability)
    .def("get_contact_force_covariance", &ContactEstimator::getContactForceCovariance)
    .def("get_contact_surface_normal", &ContactEstimator::getContactSurfaceNormal)
    .def("set_contact_surface_normal", &ContactEstimator::setContactSurfaceNormal,
          py::arg("contact_surface_normal"))
     .def("__str__", [](const ContactEstimator& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace legged_state_estimator