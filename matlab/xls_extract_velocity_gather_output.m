function[] = xls_extract_velocity_gather_output(file_name, transtion_angle)
%  file_name = 'random_S1_11gy_11gx_15000WL_1_knick_9_output';     
warning('off','all');
solid_fraction = 0.8;%965/(34*34);
g = 10;
d = 0.001;
rho_p = 2500;
rho = rho_p * solid_fraction;
depth = 34*0.001;
properties = {'vol_flux', 'velocity_profile'};

     [data, header] = xlsread(strcat('..\matlab\',file_name,'.xlsx'), 'means_velocity_profile');
     [stds, header] = xlsread(strcat('..\matlab\',file_name,'.xlsx'), 'stds_velocity_profile');
   
  from_col = 3;
  to_col = 12;
   layers = (-9:0).';
     angle_fits = cell(1, numel(data(1:end,1))+1);
     angle_plots = cell(1, numel(data(1:end,1)));
     goodnesses = cell(1, numel(data(1:end,1))+1);
     

     furb_velocity = strcat('u*log((0.5-x)/10)');
     standa_velocity = strcat(num2str(2*rho*g),'* (',num2str(sind(theta)),' - m_d*', num2str(cosd(theta)), ...
                            ')/(b* ', num2str(3*d*sqrt(rho_p/rho*g*cosd(theta))),') * (', ...
                            num2str(depth^1.5),' - ((0.5-x)*',num2str(depth/10),').^1.5) ');

     %
  % before transition angle
      for theta=1:transition_angle
          [angle_fits{theta}, goodnesses{theta}] = fit(layers, data( theta , from_col:to_col).', 'exp1');
         
      end
% after transtion angle
      for theta=1:transition_angle
          [angle_fits{theta}, goodnesses{theta}] = fit(layers, data( theta , from_col:to_col).', furb_velocity);
         
      end
      
      
      
      
      
%where x is depth inside the flow measured from the top surface
% v = 0 for y = h



      for theta=transition_angle:numel(data(1:end,1))
         % b0 = sqrt(h)*rho*g*(sind(theta) - num2str(st2_f.m_d)*cosd(theta))/ ... 
            %            (num2str(st2_f.b)*d*sqrt(rho_p*rho*g*cosd(theta)));
   
            
       %   standa_full = (2*rho*g)* (sind(theta) - st2_f.m_d*cosd(theta)) ...
          %                 /(st2_f.b* 3*d*sqrt(rho_p/rho*g*cosd(theta))) * (h^1.5) - ( h-x).^1.5;
          [angle_fits{theta+1}, goodnesses{theta+1}] = fit(layers, data( theta , from_col:to_col).', standa_velocity);
    %      angle_plot{theta}=plot(angle_fits{theta+1},layers, data( theta , from_col:to_col).');
      %    set(angle_plot{theta}(2), 'color', colors(theta,1:3));
      end
text(0.00002,  -4, r2s);
      legend off;
%view([90,-90]);
  ylim([0 2e-5])
  xlabel('depth');
  ylabel('velocity');
  
leg_name(2:2:22) = data(1:end,13);
leg_name(1:2:22)=atand(data(1:end,1));

   lh = legend(arrayfun(@num2str, leg_name, 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);     
        saveas(gcf, strcat('velocity ', file_name));


end