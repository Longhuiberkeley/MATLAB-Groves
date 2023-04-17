function  [est_C_b_e_new,est_v_eb_e_new,est_r_eb_e_new,est_IMU_bias_new,...
            est_clock_new,P_matrix_new]= NonHolocorrection(est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
            est_IMU_bias_old,est_clock_old,P_matrix_old, Vehicle_config)


velocity = est_C_b_e_old' * est_v_eb_e_old; 


reduction_matrix = [0,1,0;
    0,0,1]; 

H_nonholo = [zeros(2,3), -reduction_matrix* est_C_b_e_old', zeros(2,11)]; 

R_matrix = diag([Vehicle_config.latnoise, Vehicle_config.vertnoise]); 

K_matrix = P_matrix_old * H_nonholo' / (H_nonholo *...
    P_matrix_old * H_nonholo' + R_matrix);


delta_z = - velocity(2:3); 

x_est_new =  K_matrix * delta_z;

P_matrix_new = (eye(17) - K_matrix * H_nonholo) * P_matrix_old;

% Correct attitude, velocity, and position using (14.7-9)
est_C_b_e_new = (eye(3) - Skew_symmetric(x_est_new(1:3))) * est_C_b_e_old;
est_v_eb_e_new = est_v_eb_e_old - x_est_new(4:6);
est_r_eb_e_new = est_r_eb_e_old - x_est_new(7:9);

% Update IMU bias and GNSS receiver clock estimates
est_IMU_bias_new = est_IMU_bias_old + x_est_new(10:15);
est_clock_new = est_clock_old +  x_est_new(16:17)';


end