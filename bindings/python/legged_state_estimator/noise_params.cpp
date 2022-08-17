#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "legged_state_estimator/inekf/noise_params.hpp"


namespace legged_state_estimator {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(noise_params, m) {
  py::class_<NoiseParams>(m, "NoiseParams")
    .def(py::init<>())
    .def_property("gyro_cov", &NoiseParams::getGyroscopeCov,
                   static_cast<void (NoiseParams::*)(const Eigen::Matrix3d&)>(&NoiseParams::setGyroscopeNoise))
    .def_property("accel_cov", &NoiseParams::getAccelerometerCov,
                   static_cast<void (NoiseParams::*)(const Eigen::Matrix3d&)>(&NoiseParams::setAccelerometerNoise))
    .def_property("gyro_bias_cov", &NoiseParams::getGyroscopeBiasCov,
                   static_cast<void (NoiseParams::*)(const Eigen::Matrix3d&)>(&NoiseParams::setGyroscopeBiasNoise))
    .def_property("accel_bias_cov", &NoiseParams::getAccelerometerBiasCov,
                   static_cast<void (NoiseParams::*)(const Eigen::Matrix3d&)>(&NoiseParams::setAccelerometerBiasNoise))
    .def_property("contact_cov", &NoiseParams::getContactCov,
                   static_cast<void (NoiseParams::*)(const Eigen::Matrix3d&)>(&NoiseParams::setContactNoise));
}

} // namespace python
} // namespace legged_state_estimator