#include "inekf/state_estimator_settings.hpp"


namespace inekf {

StateEstimatorSettings StateEstimatorSettings::UnitreeA1(
    const std::string& path_to_urdf, const double dt) {
  StateEstimatorSettings settings;
  settings.path_to_urdf = path_to_urdf;
  settings.imu_frame = 46;
  settings.contact_frames = {14, 24, 34, 44}; // LF, RF, LH, RH

  settings.contact_estimator_settings.beta0 = {-20.0, -20.0, -20.0, -20.0};
  settings.contact_estimator_settings.beta1 = {0.7, 0.7, 0.7, 0.7};
  settings.contact_estimator_settings.contact_force_cov_alpha = 100.0;
  settings.contact_estimator_settings.contact_prob_threshold = 0.5;

  settings.noise_params.setGyroscopeNoise(0.01);
  settings.noise_params.setAccelerometerNoise(0.1);
  settings.noise_params.setGyroscopeBiasNoise(0.00001);
  settings.noise_params.setAccelerometerBiasNoise(0.0001);
  settings.noise_params.setContactNoise(0.1);

  settings.dynamic_contact_estimation = false;

  settings.contact_position_noise = 0.01;
  settings.contact_rotation_noise = 0.01;

  settings.dt = dt;

  settings.lpf_gyro_accel_cutoff = 250;
  settings.lpf_lin_accel_cutoff  = 250;
  settings.lpf_dqJ_cutoff        = 10;
  settings.lpf_ddqJ_cutoff       = 5;
  settings.lpf_tauJ_cutoff       = 10;

  return settings;
}

} // namespace inekf
