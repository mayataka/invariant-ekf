#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "legged_state_estimator/state_estimator_settings.hpp"


namespace legged_state_estimator {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(state_estimator_settings, m) {
  py::class_<LeggedStateEstimatorSettings>(m, "LeggedStateEstimatorSettings")
    .def(py::init<>())
    .def_static("UnitreeA1", &LeggedStateEstimatorSettings::UnitreeA1,
                 py::arg("path_to_urdf"), py::arg("dt"))
    .def_readwrite("path_to_urdf", &LeggedStateEstimatorSettings::path_to_urdf)
    .def_readwrite("contact_frames", &LeggedStateEstimatorSettings::contact_frames)
    .def_readwrite("contact_estimator", &LeggedStateEstimatorSettings::contact_estimator)
    .def_readwrite("slip_estimator_settings", &LeggedStateEstimatorSettings::slip_estimator_settings)
    .def_readwrite("noise_params", &LeggedStateEstimatorSettings::noise_params)
    .def_readwrite("dynamic_contact_estimation", &LeggedStateEstimatorSettings::dynamic_contact_estimation)
    .def_readwrite("contact_position_noise", &LeggedStateEstimatorSettings::contact_position_noise)
    .def_readwrite("contact_rotation_noise", &LeggedStateEstimatorSettings::contact_rotation_noise)
    .def_readwrite("dt", &LeggedStateEstimatorSettings::dt)
    .def_readwrite("lpf_gyro_accel_cutoff", &LeggedStateEstimatorSettings::lpf_gyro_accel_cutoff)
    .def_readwrite("lpf_lin_accel_cutoff", &LeggedStateEstimatorSettings::lpf_lin_accel_cutoff)
    .def_readwrite("lpf_dqJ_cutoff", &LeggedStateEstimatorSettings::lpf_dqJ_cutoff)
    .def_readwrite("lpf_ddqJ_cutoff", &LeggedStateEstimatorSettings::lpf_ddqJ_cutoff)
    .def_readwrite("lpf_tauJ_cutoff", &LeggedStateEstimatorSettings::lpf_dqJ_cutoff);
}

} // namespace python
} // namespace legged_state_estimator