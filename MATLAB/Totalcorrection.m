function  [est_C_b_e_new,est_v_eb_e_new,est_r_eb_e_new,est_IMU_bias_new,...
            est_clock_new,P_matrix_new]= Totalcorrection(GNSS_measurements,...
            no_meas,tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
            est_IMU_bias_old,est_clock_old,P_matrix_old, TC_KF_config, Vehicle_config)


% Constants (sone of these could be changed to inputs at a later date)
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]); % Skew symmetric matrix of Earth rate
c = 299792458; % Speed of light in m/s

% Propagate state estimates using (3.14) noting that only the clock
% states are non-zero due to closed-loop correction.
x_est_propagated = zeros(17,1); 
% x_est_propagated(1:15,1) = 0;
x_est_propagated(16,1) = est_clock_old(1) + est_clock_old(2) * tor_s;
x_est_propagated(17,1) = est_clock_old(2);


% MEASUREMENT UPDATE PHASE
u_as_e_T = zeros(no_meas,3);
pred_meas = zeros(no_meas,2);



% Loop measurements
for j = 1:no_meas

    % Predict approx range 
    delta_r = GNSS_measurements(j,3:5)' - est_r_eb_e_old;
    approx_range = sqrt(delta_r' * delta_r);

    % Calculate frame rotation during signal transit time using (8.36)
    C_e_I = [1, omega_ie * approx_range / c, 0;...
             -omega_ie * approx_range / c, 1, 0;...
             0, 0, 1];

    % Predict pseudo-range using (9.165)
    delta_r = C_e_I *  GNSS_measurements(j,3:5)' - est_r_eb_e_old;
    range = sqrt(delta_r' * delta_r);
    pred_meas(j,1) = range + x_est_propagated(16);
        
    % Predict line of sight
    u_as_e_T(j,1:3) = delta_r' / range;
        
    % Predict pseudo-range rate using (9.165)
    range_rate = u_as_e_T(j,1:3) * (C_e_I * (GNSS_measurements(j,6:8)' +...
        Omega_ie * GNSS_measurements(j,3:5)') - (est_v_eb_e_old +...
        Omega_ie * est_r_eb_e_old));        
    pred_meas(j,2) = range_rate + x_est_propagated(17);

end % for j


% Set-up measurement matrix using (14.126)
H = zeros((2 * no_meas),17);
H(1:no_meas,7:9) = u_as_e_T(1:no_meas,1:3);
H(1:no_meas,16) = ones(no_meas,1);
H((no_meas + 1):(2 * no_meas),4:6) = u_as_e_T(1:no_meas,1:3);
H((no_meas + 1):(2 * no_meas),17) = ones(no_meas,1);

velocity = est_C_b_e_old' * est_v_eb_e_old; 

reduction_matrix = [0,1,0;
    0,0,1]; 

H_nonholo = [zeros(2,3), -reduction_matrix* est_C_b_e_old', zeros(2,11)]; 

H_matrix = [H; H_nonholo]; 


% Set-up measurement noise covariance matrix assuming all measurements
% are independent and have equal variance for a given measurement type.
R_matrix(1:no_meas,1:no_meas) = eye(no_meas) *...
    TC_KF_config.pseudo_range_SD^2;
R_matrix(1:no_meas,(no_meas + 1):(2 * no_meas)) =...
    zeros(no_meas);
R_matrix((no_meas + 1):(2 * no_meas),1:no_meas) =...
    zeros(no_meas);
R_matrix((no_meas + 1):(2 * no_meas),(no_meas + 1):(2 * no_meas)) =...
    eye(no_meas) * TC_KF_config.range_rate_SD^2;

R_matrix(end+1, end+1) = Vehicle_config.latnoise;
R_matrix(end+1, end+1) = Vehicle_config.vertnoise;

% Calculate Kalman gain using (3.21)

% P_matrix_old * H_matrix' 
% [17* 17] [17* 20]
% /([20*17] [17*17]*[17*20] + R) 

K_matrix = P_matrix_old * H_matrix' / (H_matrix *...
    P_matrix_old * H_matrix' + R_matrix);



% Formulate measurement innovations using (14.119)
delta_z(1:no_meas,1) = GNSS_measurements(1:no_meas,1) -...
    pred_meas(1:no_meas,1);
delta_z((no_meas + 1):(2 * no_meas),1) = GNSS_measurements(1:no_meas,2) -...
    pred_meas(1:no_meas,2);

delta_z = [delta_z; -velocity(2:3)]; 

% Update state estimates using (3.24)
x_est_new = x_est_propagated + K_matrix * delta_z;

% Update state estimation error covariance matrix using (3.25)
P_matrix_new = (eye(17) - K_matrix * H_matrix) * P_matrix_old;

% Correct attitude, velocity, and position using (14.7-9)
est_C_b_e_new = (eye(3) - Skew_symmetric(x_est_new(1:3))) * est_C_b_e_old;
est_v_eb_e_new = est_v_eb_e_old - x_est_new(4:6);
est_r_eb_e_new = est_r_eb_e_old - x_est_new(7:9);

% Update IMU bias and GNSS receiver clock estimates
est_IMU_bias_new = est_IMU_bias_old + x_est_new(10:15);
est_clock_new = x_est_new(16:17)';



end