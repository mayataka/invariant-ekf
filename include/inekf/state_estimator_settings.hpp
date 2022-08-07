#ifndef INEKF_STATE_ESTIMATOR_SETTINGS_HPP_
#define INEKF_STATE_ESTIMATOR_SETTINGS_HPP_

#include <string>
#include <vector>

#include "inekf/macros.hpp"
#include "inekf/noise_params.hpp"
#include "inekf/contact_estimator.hpp"
#include "inekf/slip_estimator.hpp"


namespace inekf {

struct StateEstimatorSettings {
public:
  /// 
  /// @brief Path to the URDF file.
  ///
  std::string path_to_urdf;

  /// 
  /// @brief id of the IMU frame.
  ///
  int imu_frame;

  /// 
  /// @brief ids of the contact frames.
  ///
  std::vector<int> contact_frames;

  /// 
  /// @brief Contact estimator settings. 
  ///
  ContactEstimatorSettings contact_estimator;

  /// 
  /// @brief Slip estimator settings. 
  ///
  SlipEstimatorSettings slip_estimator_settings;

  /// 
  /// @brief Noise parameters (covariances) of InEKF. 
  ///
  NoiseParams noise_params;

  /// 
  /// @brief Use dynamics in contact estimation. If false, equilibrium is 
  /// used for contact estimation. Default is false.
  ///
  bool dynamic_contact_estimation = false;

  /// 
  /// @brief Noise (covariance) on contact position. (Possibly is not used in 
  /// InEKF. Contact covariance in noise_params are more important).
  ///
  double contact_position_noise;

  /// 
  /// @brief Noise (covariance) on contact rotation. Only used with surface 
  /// contacts.
  ///
  double contact_rotation_noise;

  /// 
  /// @brief Time step of estimation. 
  ///
  double dt;

  /// 
  /// @brief Cutoff frequency of LPF for gyro sensor. 
  ///
  double lpf_gyro_cutoff;

  /// 
  /// @brief Cutoff frequency of LPF for acceleration of gyro sensor. 
  ///
  double lpf_gyro_accel_cutoff;

  /// 
  /// @brief Cutoff frequency of LPF for linear acceleration measurement from IMU. 
  ///
  double lpf_lin_accel_cutoff;

  /// 
  /// @brief Cutoff frequency of LPF for joint velocities. 
  ///
  double lpf_dqJ_cutoff;

  /// 
  /// @brief Cutoff frequency of LPF for joint accelerations. 
  ///
  double lpf_ddqJ_cutoff;

  /// 
  /// @brief Cutoff frequency of LPF for joint torques. 
  ///
  double lpf_tauJ_cutoff;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  static StateEstimatorSettings UnitreeA1(const std::string& path_to_urdf, 
                                          const double dt);

};

} // namespace inekf

#endif // INEKF_STATE_ESTIMATOR_SETTINGS_HPP_