#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "legged_state_estimator/robot_model.hpp"


namespace legged_state_estimator {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(robot_model, m) {
  py::enum_<pinocchio::ReferenceFrame>(m, "ReferenceFrame")
    .value("WORLD", pinocchio::ReferenceFrame::WORLD)
    .value("LOCAL", pinocchio::ReferenceFrame::LOCAL)
    .value("LOCAL_WORLD_ALIGNED", pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED)
    .export_values();

  py::class_<RobotModel>(m, "RobotModel")
    .def(py::init<const std::string&, const int, const std::vector<int>&>(),
          py::arg("path_to_urdf"), py::arg("imu_frame"), py::arg("contact_frames"))
    .def(py::init<const std::string&, const std::string&, const std::vector<std::string>&>(),
          py::arg("path_to_urdf"), py::arg("imu_frame"), py::arg("contact_frames"))
    .def(py::init<>())
    .def("update_leg_kinematics", static_cast<void (RobotModel::*)(const Eigen::VectorXd&, 
                                                                   const pinocchio::ReferenceFrame)>(&RobotModel::updateLegKinematics),
          py::arg("qJ"), py::arg("rf")=pinocchio::LOCAL_WORLD_ALIGNED)
    .def("update_leg_kinematics", static_cast<void (RobotModel::*)(const Eigen::VectorXd&, const Eigen::VectorXd&, 
                                                                   const pinocchio::ReferenceFrame)>(&RobotModel::updateLegKinematics),
          py::arg("qJ"), py::arg("dqJ"), py::arg("rf")=pinocchio::LOCAL_WORLD_ALIGNED)
    .def("update_kinematics", static_cast<void (RobotModel::*)(const Eigen::Vector3d&, 
                                                               const Eigen::Vector4d&,  
                                                               const Eigen::VectorXd&)>(&RobotModel::updateKinematics),
          py::arg("base_pos"), py::arg("base_quat"), py::arg("qJ"))
    .def("update_kinematics", static_cast<void (RobotModel::*)(const Eigen::Vector3d&, 
                                                               const Eigen::Vector4d&, 
                                                               const Eigen::Vector3d&, 
                                                               const Eigen::Vector3d&, 
                                                               const Eigen::VectorXd&, 
                                                               const Eigen::VectorXd&)>(&RobotModel::updateKinematics),
          py::arg("base_pos"), py::arg("base_quat"), 
          py::arg("base_linear_vel"), py::arg("base_angular_vel"), 
          py::arg("qJ"), py::arg("dqJ"))
    .def("update_leg_dynamics", &RobotModel::updateLegDynamics,
          py::arg("qJ"), py::arg("dqJ"))
    .def("update_dynamics", &RobotModel::updateDynamics,
          py::arg("base_pos"), py::arg("base_quat"), 
          py::arg("base_linear_vel"), py::arg("base_angular_vel"), 
          py::arg("base_linear_acc"), py::arg("base_angular_acc"), 
          py::arg("qJ"), py::arg("dqJ"), py::arg("ddqJ"))
    .def("get_base_position", &RobotModel::getBasePosition)
    .def("get_base_rotation", &RobotModel::getBaseRotation)
    .def("get_contact_position", &RobotModel::getContactPosition,
          py::arg("contact_id"))
    .def("get_contact_rotation", &RobotModel::getContactRotation,
          py::arg("contact_id"))
    .def("get_contact_jacobian", &RobotModel::getContactJacobian,
          py::arg("contact_id"))
    .def("get_joint_contact_jacobian", &RobotModel::getJointContactJacobian,
          py::arg("contact_id"))
    .def("get_inverse_dynamics", &RobotModel::getInverseDynamics)
    .def("get_joint_inverse_dynamics", &RobotModel::getJointInverseDynamics);
}

} // namespace python
} // namespace legged_state_estimator