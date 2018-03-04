function[] = plot_from_output_gg(dir, file_name, transition_angle, num_of_layers)
%       
warning('off','all');
solid_fraction = 965/(38*34);%0.8
g = 10;
d = 1;
rho_p = 2500;
rho = rho_p * solid_fraction;
h = 38;
dt1 = 0.07158707624543059;
dt2 = 7.1570724900297727E-002;
dt = dt1;
set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
colors = get(0, 'DefaultAxesColorOrder');
from_col = 3;
to_col = num_of_layers + 2;
total_col = 13;%to_col + 1;
properties = {'vol_flux', 'velocity_profile'};
%file_name = 'random_S1_11gy_11gx_15000WL_1_knick_9_output';
     [flux_data, header] = xlsread(strcat(dir, file_name,'.xlsx'), 'means_vol_flux');
     [flux_stds, header] = xlsread(strcat(dir, file_name,'.xlsx'), 'stds_vol_flux');
     [vel_data, header] = xlsread(strcat(dir, file_name,'.xlsx'), 'means_velocity_profile');
     [vel_stds, header] = xlsread(strcat(dir, file_name,'.xlsx'), 'stds_velocity_profile');
     output_fig = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', file_name, '_velocity_plots');
       fig_H  = figure('name', file_name);
         title ('1 with solfrac');
       legend off;
       % u(z=Hs) = 0
     Hs = h+h/(2*(num_of_layers));
     num_of_angles = numel(vel_data(1:end,1));
    standa = strcat(num2str(2 * sqrt(solid_fraction)),'*(cosd(x)^0.5) *(tand(x) - m_d) * ', num2str(3/5*(Hs)^(5/2)  / (3 *d)),'/b');
    furbish = strcat(num2str(Hs),'*a*tand(x)/((1-tand(x)/b)*cosd(x))');
    furbish_2 = strcat(num2str(Hs/rho/g/d^4),...
        '*a/alpha_cos_phi*tand(x)/((1-tand(x)*b/alpha_cos_phi/',num2str(d^2),')*cosd(x))');
        %axis([0,0,tand(16),10]) ;
    
        %%%%  start flux  %%%%
        subplot(2,1,1);
        hold on;        
  
    [furb1_f, furb1_goodness] = fit(flux_data(1:transition_angle,1), flux_data(1:transition_angle, total_col), furbish,  'lower', [0, 0.1], 'upper', [10, 0.5]);
    while furb1_goodness.rsquare < 0.95
        display(strcat('bad furbish fit-  R^2 = ', num2str(furb1_goodness.rsquare)));
        furb1_f
        [furb1_f, furb1_goodness] = fit(flux_data(1:transition_angle,1), flux_data(1:transition_angle, total_col), furbish,  'lower', [0, 0.1], 'upper', [10, 0.5]);
    end
    furb1_plot = plot(furb1_f, flux_data(1:transition_angle,1), flux_data(1:transition_angle, total_col));

    set(furb1_plot(2), 'color', 'k');

    [st2_f, st2_goodness] = fit(flux_data(transition_angle:end,1), flux_data(transition_angle:end, total_col),standa);
    st2_plot = plot(st2_f, flux_data(transition_angle:end,1), flux_data(transition_angle:end, total_col));
    set(st2_plot(2), 'color', 'b');
     e = errorbar(flux_data(1:end,1),flux_data(1:end, total_col),flux_stds(1:end, total_col), '.'); 
     e.LineStyle = 'none';
    title('with L & solfrac');
    xlabel('slope angle');  
    ylabel('average total flux');
    legend('location', 'northwest');
        
    text((flux_data(1,1)+flux_data(floor(num_of_angles/2),1))/4,flux_data(end, total_col),...
          {strcat('Furbish (k): K=',num2str(furb1_f.a), ' b=', num2str(furb1_f.b), ' Sc (atand(b))=', num2str(atand(furb1_f.b)),' rsquare=', num2str(furb1_goodness.rsquare)),...                 
                    strcat('Standa (b): m_d=',num2str(st2_f.m_d), ' beta=', num2str(st2_f.b),' rsquare=', num2str(st2_goodness.rsquare))},...
                    ...%strcat('Furbish (r): K=',num2str(furb_f.a), ' Sc=', num2str(atand(furb_f.b)),' rsquare=', num2str(furb_goodness.rsquare))},...
          'HorizontalAlignment','left',...
          'VerticalAlignment','top',...
          'BackgroundColor',[1 1 1],...
          'FontSize',10);
     
    legend('off');
    
  
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % depth by mid sub layers
depth =  -(Hs : -h/(num_of_layers) : h/(2*(num_of_layers))).';
% bottom of layers depth = -(h : -h/(num_of_layers) : h/(num_of_layers)).';
h0=-Hs/exp(1); % == 0.3627*depth(7)+(1-0.3627)*depth(8)  

% figure1=  figure('name', strcat(sub_name, '_vel_fit'));
   %   axes1 = axes('Parent',figure1);
%view(axes1,[90 -90]);
%hold(axes1,'all');

  sub_ax = subplot(2,1,2);
 
  %view(sub_ax,[90 -90]);
        hold on;   
     
        
        plot(depth, ([zeros(num_of_angles,1) vel_data(1:end, from_col:to_col)]).', '.');
    vels = [zeros(num_of_angles,1) vel_data(1:end, from_col:to_col)];
 
    
    % + 1 for the transition angle (fit both furbish and standa)
    curve_from_data=zeros(num_of_angles+1, num_of_layers+1);
    curve_from_fit=zeros(num_of_angles+1, num_of_layers+1);
    u0_from_fit=zeros(num_of_angles+1,1);
    u0_from_data=zeros(num_of_angles+1,1);
    h0s=zeros(num_of_angles+1,1);
    Hss_idx=ones(num_of_angles+1,1);
    
     angle_fits = cell(1, num_of_angles+1);
     angle_plots = cell(1, num_of_angles);
     goodnesses = cell(1, num_of_angles+1);
     angle_plots_out_a = cell(1, num_of_angles);
     angle_plots_out_b = cell(1, num_of_angles);
     goodnesses_out = cell(1, num_of_angles+1);
     
     angle_plots_out_f1 = cell(1, num_of_angles);
     angle_plots_out_f2 = cell(1, num_of_angles);
     angle_plots_out_f3 = cell(1, num_of_angles);
    
     % before transition angle
    r2s='furbish::';
      for theta=1:transition_angle
         
         tan_th = tand(vel_data(theta,1));
         cos_th = cosd(vel_data(theta,1));
         u0_from_data(theta)=0.3627*vels(theta,7)+(1-0.3627)*vels(theta,8); % assuming 10 sub layers
         
         
         u0_from_fit(theta) = furb1_f.a/(dt*50000)/solid_fraction/0.7854*tan_th/((1-tan_th/furb1_f.b)*cos_th);
         
         % find the depth corresponding to u0 from fit
         fraction = 100;
         for layer = 2:num_of_layers+1
             if vels(theta, layer) >= u0_from_fit(theta)
                 % this theta's h0 is between this layer and the previous one
                 syms x;
                 fraction = eval(solve(x*vels(theta, layer-1)+(1-x)*vels(theta, layer) == u0_from_fit(theta)));
                 h0s(theta) = (fraction*depth(layer-1) + (1 - fraction)*depth(layer));
                 
                 % find the depth index for this angle 
                 %(where vel = 0 according to curr angle h0)
                 while depth(Hss_idx(theta)) < h0s(theta)*exp(1)
                     Hss_idx(theta) = Hss_idx(theta) + 1;
                 end
                 
                 break;
             end
         end    
         
          % with top layer (no top is num_of_layers+1)
          curve_from_data(theta, 1:end) = u0_from_data(theta)*(1-log(depth/h0));
          res = sum((vels(theta, 1:end) -  curve_from_data(theta, 1:end)).^2);
          tot = sum((vels(theta, 1:end) -  mean(vels(theta, 1:end))).^2);
          rsq = 1 - res/tot;
          display(strcat('from data theta =  ', num2str(vel_data(theta,1)), ' R^2 = ', num2str(rsq)));
          
          %curve_from_fit(theta, 1:end) =u0_from_fit(theta)*(1-log(depth/h0));
          
          if fraction < 100
              curve_from_fit(theta, Hss_idx(theta):end) =u0_from_fit(theta)*(1-log(depth(Hss_idx(theta):end)./h0s(theta)));

              res = sum((vels(theta, 1:end) -  curve_from_fit(theta, 1:end)).^2);
              rsq = 1 - res/tot;
              display(strcat('from furbish params theta =  ', num2str(vel_data(theta,1)), ' R^2 = ', num2str(rsq)));
          else
              display(strcat('no h0 for theta =  ', num2str(vel_data(theta,1)), 'u0_from_fit ',  ...
                   num2str(u0_from_fit(theta)), ' and top velocity  ', num2str(vels(theta, end))));
          end
          
          
          
            %r2s=strcat(' ', r2s, num2str(theta), ' : ', num2str(goodnesses{theta}.rsquare));
          %set(angle_plots{theta}(2), 'color', colors(theta,1:3));
          
      end           
     % plot(depth(2:end), (curve_from_fit(1:transition_angle,2:end)).');
       plot(depth, (curve_from_data(1:transition_angle,1:end)).');      
      
% after transition angle                        
      for theta=transition_angle:num_of_angles

          sin_th =sind(vel_data(theta,1));
          cos_th = cosd(vel_data(theta,1));
          standa_for_fit = strcat(num2str(2),'* (',num2str(sin_th),' - m_d*', num2str(cos_th), ...
                        ')/(b* ', num2str(3*d*sqrt(rho_p/rho*cos_th)),') * (',num2str(Hs^1.5),' - (-x).^1.5) ');
         standa_from_fit = 2*(sin_th - st2_f.m_d*cos_th)/(st2_f.b*3*d*sqrt(rho_p/rho*cos_th));
        curve_from_fit(theta+1, 1:end) =standa_from_fit*(Hs^1.5 - (-depth).^1.5);
        res = sum((vels(theta, 1:end) -  curve_from_fit(theta+1, 1:end)).^2);
      tot = sum((vels(theta, 1:end) -  mean(vels(theta, 1:end))).^2);
      rsq = 1 - res/tot;
      display(strcat('from standa params theta =  ', num2str(vel_data(theta,1)), ' R^2 = ', num2str(rsq)));
%where x is depth inside the flow measured from the top surface
% v = 0 for y = Hs
      
      end
      text(0.00002,  -4, r2s);

     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
  
% names
%leg_name(2:2:22) = vel_data(1:end,total_col);
%leg_name(1:2:22)=atand(vel_data(1:end,1));
%plot(depth, curve_from_data.');
%plot(depth, curve_from_fit.');
plot(depth, ([zeros(num_of_angles,1) vel_data(1:end, from_col:to_col)]).', '.');
      view(sub_ax,[90 -90]);
      ylabel('average velocity');  
    xlabel('depth');
 %  lh = legend(arrayfun(@num2str, leg_name, 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);     
 lh = legend(arrayfun(@num2str, vel_data(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);     
      saveas(gcf, strcat(strrep(output_fig, '.','-'), '_basic'));
     
        
        output_file = strcat('output/',file_name,'_velocity_fit.xlsx');
        
    % xlswrite(output_file, vel_data(1:end,1).', 'velocity_fits', 'B1');
% xlswrite(output_file, {'furbish AKf/c'}, 'velocity_fits', 'A2');
% xlswrite(output_file, {'furbish alpha_gamma'}, 'velocity_fits', 'A3');
% xlswrite(output_file, {'furbish alpha*cos(phi)'}, 'velocity_fits', 'A4');
% xlswrite(output_file, {'furbish R^2'}, 'velocity_fits', 'A5');
% xlswrite(output_file, {'Standa m_d'}, 'velocity_fits', 'A6');
% xlswrite(output_file, {'Standa beta'}, 'velocity_fits', 'A7');
% xlswrite(output_file, {'Standa R^2'}, 'velocity_fits', 'A8');
% xlswrite(output_file, {angle_plots_out_f1{(1:transition_angle)}}, 'velocity_fits', 'B2');
% xlswrite(output_file, {angle_plots_out_f2{(1:transition_angle)}}, 'velocity_fits', 'B3');
% xlswrite(output_file, {angle_plots_out_f3{(1:transition_angle)}}, 'velocity_fits', 'B4');
% xlswrite(output_file, {goodnesses_out{(1:transition_angle)}}, 'velocity_fits', 'B5');
% xlswrite(output_file, {angle_plots_out_a{(transition_angle+1:num_of_angles+1)}}, 'velocity_fits', strcat(char(65+transition_angle),'6'));
% xlswrite(output_file, {angle_plots_out_b{(transition_angle+1:num_of_angles+1)}}, 'velocity_fits', strcat(char(65+transition_angle),'7'));
% xlswrite(output_file, {goodnesses_out{(transition_angle+1:num_of_angles+1)}}, 'velocity_fits', strcat(char(65+transition_angle),'8'));

%      
% xlswrite(output_file, vel_data(1:end,1).', 'velocity_fits', 'B1');
% xlswrite(output_file, {'u0_from_data'}, 'velocity_fits', 'A2');
% xlswrite(output_file, {'u0_from_fit'}, 'velocity_fits', 'A3');
% xlswrite(output_file, {'furbish K'}, 'velocity_fits', 'A3');
% xlswrite(output_file, {'furbish Sc'}, 'velocity_fits', 'A4');
% xlswrite(output_file, {'furbish R^2'}, 'velocity_fits', 'A5');
% xlswrite(output_file, {'Standa m_d'}, 'velocity_fits', 'A6');
% xlswrite(output_file, {'Standa beta'}, 'velocity_fits', 'A7');
% xlswrite(output_file, {'Standa R^2'}, 'velocity_fits', 'A8');
% xlswrite(output_file, u0_from_data(1:transition_angle), 'velocity_fits', 'B2');
% xlswrite(output_file, u0_from_fit(1:transition_angle), 'velocity_fits', 'B3');
% xlswrite(output_file, {angle_plots_out_a{(1:transition_angle)}}, 'velocity_fits', 'B4');
% xlswrite(output_file, {angle_plots_out_b{(1:transition_angle)}}, 'velocity_fits', 'B5');
% xlswrite(output_file, {goodnesses_out{(1:transition_angle)}}, 'velocity_fits', 'B6');
% xlswrite(output_file, {angle_plots_out_a{(transition_angle+1:num_of_angles+1)}}, 'velocity_fits', strcat(char(65+transition_angle),'7'));
% xlswrite(output_file, {angle_plots_out_b{(transition_angle+1:num_of_angles+1)}}, 'velocity_fits', strcat(char(65+transition_angle),'8'));
% xlswrite(output_file, {goodnesses_out{(transition_angle+1:num_of_angles+1)}}, 'velocity_fits', strcat(char(65+transition_angle),'9'));

end