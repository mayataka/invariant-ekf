#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "inekf/state_estimator.hpp"


namespace inekf {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(state_estimator, m) {
  py::class_<StateEstimator>(m, "StateEstimator")
    .def(py::init<const StateEstimatorSettings&>(),
          py::arg("state_estimator_settings"))
    .def(py::init<>())
    .def("init", static_cast<void (StateEstimator::*)(const Eigen::Vector3d&, 
                                                      const Eigen::Vector4d&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&)>(&StateEstimator::init),
          py::arg("base_pos"), py::arg("base_quat"), 
          py::arg("base_lin_vel_world")=Eigen::Vector3d::Zero(), 
          py::arg("imu_gyro_bias")=Eigen::Vector3d::Zero(), 
          py::arg("imu_lin_accel_bias")=Eigen::Vector3d::Zero())
    .def("init", static_cast<void (StateEstimator::*)(const Eigen::Vector3d&, 
                                                      const Eigen::Vector4d&,
                                                      const Eigen::VectorXd&,
                                                      const std::vector<double>&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&,
                                                      const Eigen::Vector3d&)>(&StateEstimator::init),
          py::arg("base_pos"), py::arg("base_quat"), py::arg("qJ"), 
          py::arg("ground_height")=std::vector<double>({0., 0., 0., 0.}), 
          py::arg("base_lin_vel_world")=Eigen::Vector3d::Zero(), 
          py::arg("imu_gyro_bias")=Eigen::Vector3d::Zero(), 
          py::arg("imu_lin_accel_bias")=Eigen::Vector3d::Zero())
    .def("update", &StateEstimator::update,
          py::arg("imu_gyro_raw"), py::arg("imu_lin_accel_raw"), 
          py::arg("qJ"), py::arg("dqJ"), py::arg("tauJ"))
    .def("get_contact_estimator", &StateEstimator::getContactEstimator)
    .def("get_slip_estimator", &StateEstimator::getSlipEstimator)
    .def("reset_contact_surface_normal_estimate", &StateEstimator::resetContactSurfaceNormalEstimate,
          py::arg("contact_surface_normal"))
    .def("reset_friction_coefficient_estimate", &StateEstimator::resetFrictionCoefficientEstimate,
          py::arg("friction_coefficient"))
    .def_property_readonly("base_position_estimate", &StateEstimator::getBasePositionEstimate)
    .def_property_readonly("base_rotation_estimate", &StateEstimator::getBaseRotationEstimate)
    .def_property_readonly("base_quaternion_estimate", &StateEstimator::getBaseQuaternionEstimate)
    .def_property_readonly("base_linear_velocity_estimate_world", &StateEstimator::getBaseLinearVelocityEstimateWorld)
    .def_property_readonly("base_linear_velocity_estimate_local", &StateEstimator::getBaseLinearVelocityEstimateLocal)
    .def_property_readonly("base_angular_velocity_estimate_world", &StateEstimator::getBaseAngularVelocityEstimateWorld)
    .def_property_readonly("base_angular_velocity_estimate_local", &StateEstimator::getBaseAngularVelocityEstimateLocal)
    .def_property_readonly("imu_gyro_bias_estimate", &StateEstimator::getIMUGyroBiasEstimate)
    .def_property_readonly("imu_linear_acceleration_bias_estimate", &StateEstimator::getIMULinearAccelerationBiasEstimate)
    .def_property_readonly("joint_velocity_estimate", &StateEstimator::getJointVelocityEstimate)
    .def_property_readonly("joint_acceleration_estimate", &StateEstimator::getJointAccelerationEstimate)
    .def_property_readonly("joint_torque_estimate", &StateEstimator::getJointTorqueEstimate);
}

} // namespace python
} // namespace inekf