#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "legged_state_estimator/state_estimator.hpp"


namespace legged_state_estimator {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(state_estimator, m) {
  py::class_<LeggedStateEstimator>(m, "LeggedStateEstimator")
    .def(py::init<const LeggedStateEstimatorSettings&>(),
          py::arg("state_estimator_settings"))
    .def(py::init<>())
    .def("init", static_cast<void (LeggedStateEstimator::*)(const Eigen::Vector3d&, 
                                                      const Eigen::Vector4d&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&)>(&LeggedStateEstimator::init),
          py::arg("base_pos"), py::arg("base_quat"), 
          py::arg("base_lin_vel_world")=Eigen::Vector3d::Zero(), 
          py::arg("imu_gyro_bias")=Eigen::Vector3d::Zero(), 
          py::arg("imu_lin_accel_bias")=Eigen::Vector3d::Zero())
    .def("init", static_cast<void (LeggedStateEstimator::*)(const Eigen::Vector3d&, 
                                                      const Eigen::Vector4d&,
                                                      const Eigen::VectorXd&,
                                                      const std::vector<double>&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&)>(&LeggedStateEstimator::init),
          py::arg("base_pos"), py::arg("base_quat"), py::arg("qJ"), 
          py::arg("ground_height")=std::vector<double>({0., 0., 0., 0.}), 
          py::arg("base_lin_vel_world")=Eigen::Vector3d::Zero(), 
          py::arg("imu_gyro_bias")=Eigen::Vector3d::Zero(), 
          py::arg("imu_lin_accel_bias")=Eigen::Vector3d::Zero())
    .def("update", &LeggedStateEstimator::update,
          py::arg("imu_gyro_raw"), py::arg("imu_lin_accel_raw"), 
          py::arg("qJ"), py::arg("dqJ"), py::arg("tauJ"))
    .def("get_contact_estimator", &LeggedStateEstimator::getContactEstimator)
    .def("get_slip_estimator", &LeggedStateEstimator::getSlipEstimator)
    .def("reset_contact_surface_normal_estimate", &LeggedStateEstimator::resetContactSurfaceNormalEstimate,
          py::arg("contact_surface_normal"))
    .def("reset_friction_coefficient_estimate", &LeggedStateEstimator::resetFrictionCoefficientEstimate,
          py::arg("friction_coefficient"))
    .def_property_readonly("base_position_estimate", &LeggedStateEstimator::getBasePositionEstimate)
    .def_property_readonly("base_rotation_estimate", &LeggedStateEstimator::getBaseRotationEstimate)
    .def_property_readonly("base_quaternion_estimate", &LeggedStateEstimator::getBaseQuaternionEstimate)
    .def_property_readonly("base_linear_velocity_estimate_world", &LeggedStateEstimator::getBaseLinearVelocityEstimateWorld)
    .def_property_readonly("base_linear_velocity_estimate_local", &LeggedStateEstimator::getBaseLinearVelocityEstimateLocal)
    .def_property_readonly("base_angular_velocity_estimate_world", &LeggedStateEstimator::getBaseAngularVelocityEstimateWorld)
    .def_property_readonly("base_angular_velocity_estimate_local", &LeggedStateEstimator::getBaseAngularVelocityEstimateLocal)
    .def_property_readonly("imu_gyro_bias_estimate", &LeggedStateEstimator::getIMUGyroBiasEstimate)
    .def_property_readonly("imu_linear_acceleration_bias_estimate", &LeggedStateEstimator::getIMULinearAccelerationBiasEstimate)
    .def_property_readonly("joint_velocity_estimate", &LeggedStateEstimator::getJointVelocityEstimate)
    .def_property_readonly("joint_acceleration_estimate", &LeggedStateEstimator::getJointAccelerationEstimate)
    .def_property_readonly("joint_torque_estimate", &LeggedStateEstimator::getJointTorqueEstimate);
}

} // namespace python
} // namespace legged_state_estimator