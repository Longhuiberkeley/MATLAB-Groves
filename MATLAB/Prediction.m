function [est_r_eb_e,est_v_eb_e,est_C_b_e,P_matrix_propagated] = Prediction(tor_i,est_C_b_e_old,est_r_eb_e_old,est_v_eb_e_old,...
    est_IMU_bias, P_matrix_old,meas_f_ib_b,meas_omega_ib_b,...
    est_L_b_old,TC_KF_config)

% Modified by Long Hui 

% this given, delta t, the state, imu measurement, 
% predict the next time-step state, and propagate the covariance matrix

% Inputs:
%   tor_i                 propagation interval (s)
%   est_C_b_e_old         prior estimated body to ECEF coordinate
%                         transformation matrix
%   est_v_eb_e_old        prior estimated ECEF user position (m)
%   est_r_eb_e_old        prior estimated ECEF user position (m)
%   est_IMU_bias      prior estimated IMU biases (body axes)
%   P_matrix_old          previous Kalman filter error covariance matrix
%   meas_f_ib_b           measured specific force
%   meas_omega_ib_b       meausred angular velocity
%   est_L_b_old           previous latitude solution
%   TC_KF_config
%     .gyro_noise_PSD     Gyro noise PSD (rad^2/s)
%     .accel_noise_PSD    Accelerometer noise PSD (m^2 s^-3)
%     .accel_bias_PSD     Accelerometer bias random walk PSD (m^2 s^-5)
%     .gyro_bias_PSD      Gyro bias random walk PSD (rad^2 s^-3)
%     .clock_freq_PSD     Receiver clock frequency-drift PSD (m^2/s^3)
%     .clock_phase_PSD    Receiver clock phase-drift PSD (m^2/s)
%     .pseudo_range_SD    Pseudo-range measurement noise SD (m)
%     .range_rate_SD      Pseudo-range rate measurement noise SD (m/s)
%
% Outputs:
%   est_C_b_e_new     updated estimated body to ECEF coordinate 
%                      transformation matrix
%   est_v_eb_e_new    updated estimated ECEF user velocity (m/s)
%   est_r_eb_e_new    updated estimated ECEF user position (m)
%   est_IMU_bias_new  updated estimated IMU biases



% Constants (sone of these could be changed to inputs at a later date)
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity

% Begins

% Skew symmetric matrix of Earth rate
Omega_ie = Skew_symmetric([0,0,omega_ie]);

% Correct IMU errors
meas_f_ib_b = meas_f_ib_b - est_IMU_bias(1:3);
meas_omega_ib_b = meas_omega_ib_b - est_IMU_bias(4:6);

% Update estimated navigation solution
[est_r_eb_e,est_v_eb_e,est_C_b_e] = Nav_equations_ECEF(tor_i,...
    est_r_eb_e_old,est_v_eb_e_old,est_C_b_e_old,meas_f_ib_b,...
    meas_omega_ib_b);



Phi_matrix = eye(17);
Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * tor_i;
Phi_matrix(1:3,13:15) = est_C_b_e_old * tor_i;
Phi_matrix(4:6,1:3) = -tor_i * Skew_symmetric(est_C_b_e_old * meas_f_ib_b);
Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) - 2 * Omega_ie * tor_i;
geocentric_radius = R_0 / sqrt(1 - (e * sin(est_L_b_old))^2) *...
    sqrt(cos(est_L_b_old)^2 + (1 - e^2)^2 * sin(est_L_b_old)^2); % from (2.137)
Phi_matrix(4:6,7:9) = -tor_i * 2 * Gravity_ECEF(est_r_eb_e_old) /...
    geocentric_radius * est_r_eb_e_old' / sqrt (est_r_eb_e_old' *...
    est_r_eb_e_old);
Phi_matrix(4:6,10:12) = est_C_b_e_old * tor_i;
Phi_matrix(7:9,4:6) = eye(3) * tor_i;
Phi_matrix(16,17) = tor_i;


Q_prime_matrix = zeros(17);
Q_prime_matrix(1:3,1:3) = eye(3) * TC_KF_config.gyro_noise_PSD * tor_i;
Q_prime_matrix(4:6,4:6) = eye(3) * TC_KF_config.accel_noise_PSD * tor_i;
Q_prime_matrix(10:12,10:12) = eye(3) * TC_KF_config.accel_bias_PSD * tor_i;
Q_prime_matrix(13:15,13:15) = eye(3) * TC_KF_config.gyro_bias_PSD * tor_i;
Q_prime_matrix(16,16) = TC_KF_config.clock_phase_PSD * tor_i;
Q_prime_matrix(17,17) = TC_KF_config.clock_freq_PSD * tor_i;

P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *...
    Phi_matrix' + 0.5 * Q_prime_matrix;



end