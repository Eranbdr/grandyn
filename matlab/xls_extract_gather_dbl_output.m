function[trans_angle, trans_velocity, trans_flux, K, Sc, Rr2, Kf, N, Fr2, mu_d, beta, Sr2, flux_pairs1, flux_pairs2] = xls_extract_gather_dbl_output(file_name, transition_angle)
  %transition_angle = 9;
%file_name = 'random_S1_11gy_11gx_15000WL_1_knick_9_output';
%       plot a range of columns as a single fit
warning('off','all');
solid_fraction = 0.8;%965/(34*34);
g = 10;
d = 0.001;
rho_p = 2500;
rho = rho_p * solid_fraction;
depth = 34*0.001;

% layers
layer1_col = 3;
surface_layer_col = 12;

% total
total_column = 13;

prop = 1;
properties = {'vol_flux', 'velocity_profile'};
%file_name = 'random_S1_11gy_11gx_15000WL_1_knick_9_output';
     [data, header] = xlsread(strcat('..\matlab\',file_name), strcat('means_', properties{1}));
     [stds, header] = xlsread(strcat('..\matlab\',file_name), strcat('stds_',properties{1}));
     [vel_data, header] = xlsread(strcat('..\matlab\',file_name), strcat('means_', properties{2}));

     thetas = data(1:end,1);
     flux_fit = cell(1, 4);
     flux_goodness = cell(1, 4);
     flux_plots = cell(1,4);
     
     fit_names = {'furbish all', 'furbish before', 'roering before', 'standa after'};
     furb_vel = 'u*log((0.5-x)/10)';
     furb_flux = '(a*tand(x))/((1-tand(x)/b)*cosd(tand(x)))';
     roering_flux = '(a*tand(x))/(1-(tand(x)/b)^2)';
     standa_flux = strcat(num2str(2 * rho * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_d) * ', ...
         num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *b)');
     % where x is theta
     
     trans_angle_s1 = 2 * transition_angle - 1;
     trans_angle_s2 = 2 * transition_angle;
           %fit(means(1:knick,1), means(1:knick, col),'(a*tand(x))/(1-(tand(x)/b)^2)', 'lower', [0, 0.26], 'upper', [500, 0.3]);
      [flux_fit{1}, flux_goodness{1}] = fit(thetas(1:2:end), data( 1:2:end , total_column), furb_flux);   
      
        % before transition angle
      [flux_fit{2}, flux_goodness{2}] = fit(thetas(1:2:trans_angle_s1), data( 1:2:trans_angle_s1 , total_column), furb_flux);
      
        % before transition angle
      [flux_fit{3}, flux_goodness{3}] = fit(thetas(1:2:trans_angle_s1), data( 1:2:trans_angle_s1 , total_column), roering_flux, 'lower', [0, 0.26], 'upper', [500, 0.3]);
      
      % after transtion angle
      
      [flux_fit{4}, flux_goodness{4}] = fit(thetas(trans_angle_s1:2:end), data( trans_angle_s1:2:end , total_column), standa_flux);
      
       [flux_fit{5}, flux_goodness{5}] = fit(thetas(2:2:end), data( 2:2:end , total_column), furb_flux);   
      
        % before transition angle
      [flux_fit{6}, flux_goodness{6}] = fit(thetas(2:2:trans_angle_s2), data( 2:2:trans_angle_s2 , total_column), furb_flux);
      
        % before transition angle
      [flux_fit{7}, flux_goodness{7}] = fit(thetas(2:2:trans_angle_s2), data( 2:2:trans_angle_s2 , total_column), roering_flux);
      
      % after transtion angle
      
      [flux_fit{8}, flux_goodness{8}] = fit(thetas(trans_angle_s2:2:end), data( trans_angle_s2:2:end , total_column), standa_flux);
   
  %trans_flux K Sc Rr2 Kf N Fr2 mu_d beta Sr2 flux_pairs
  trans_flux(1) = data( trans_angle_s1 , total_column);
  trans_flux(2) = data( trans_angle_s2 , total_column);
  flux_pairs1 = strcat('[', num2str(thetas(1)), ' : ', num2str(data(1,total_column), '%4.2e\n'), ']');
  flux_pairs2 = strcat('[', num2str(thetas(2)), ' : ', num2str(data(2,total_column), '%4.2e\n'), ']');
  for i = 3:2:numel(thetas)
      
      flux_pairs1 = strcat(flux_pairs1, ', [', num2str(thetas(i)), ' : ', num2str(data(i,total_column), '%4.2e\n'), ']');
      flux_pairs2 = strcat(flux_pairs2, ', [', num2str(thetas(i+1)), ' : ', num2str(data(i+1,total_column), '%4.2e\n'), ']');
  end
  K(1) = flux_fit{3}.a;
  Sc(1) = flux_fit{3}.b;
  Rr2(1) = flux_goodness{3}.rsquare;
  Kf(1) = flux_fit{2}.a;
  N(1) = flux_fit{2}.b;
  Fr2(1) = flux_goodness{2}.rsquare;
  mu_d(1) = flux_fit{4}.m_d;
  beta(1) = flux_fit{4}.b;
  Sr2(1) = flux_goodness{4}.rsquare;
  trans_velocity(1) = vel_data(trans_angle_s1 , total_column);
  trans_angle(1) = thetas(trans_angle_s1);
  
  
   K(2) = flux_fit{7}.a;
  Sc(2) = flux_fit{7}.b;
  Rr2(2) = flux_goodness{7}.rsquare;
  Kf(2) = flux_fit{6}.a;
  N(2) = flux_fit{6}.b;
  Fr2(2) = flux_goodness{6}.rsquare;
  mu_d(2) = flux_fit{8}.m_d;
  beta(2) = flux_fit{8}.b;
  Sr2(2) = flux_goodness{8}.rsquare;
  trans_velocity(2) = vel_data(trans_angle_s2 , total_column);
  trans_angle(2) = thetas(trans_angle_s2);
end