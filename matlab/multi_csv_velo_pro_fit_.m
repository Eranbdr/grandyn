function[] = multi_csv_flux_velo_pro_fit_(sub_name)

% read several CSVs and present means and vars 
% input is slope, file, SS_step
warning('off','all');
solid_fraction = 0.8;%965/(34*34);
g = 10;
d = 0.001;
rho_p = 2500;
rho = rho_p * solid_fraction;
h = 34*0.001;

columns = 13;% S 13, M 23

r_knick = 8; % from lowest slope
s_knick = 8;
set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
colors = get(0, 'DefaultAxesColorOrder');

properties = {'vol_flux', 'velocity_profile'};

[num, files] = xlsread(strcat('C:\Users\EranBD\Documents\thesis\matlab\', sub_name, '_input.csv'));
output_file = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', sub_name, '_knick_', num2str(r_knick),  '_output.xlsx');
output_fig = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', sub_name, '_knick_', num2str(r_knick));
delete output_file;
files_num = numel(files);
figure1 = figure('name', strcat(sub_name, '_knick_', num2str(r_knick)));

%%%%%  start flux    %%%%% 
%     
    [header, means_t, vars_t, stds_t] = single_csv_data(files{1}, num(1,3), num(1,4), properties{1}, 50000, columns);
    cols = numel(header);
    means= zeros(files_num, cols);
    vars = zeros(files_num, cols);
    stds = zeros(files_num, cols);
    means(1,1:end) = means_t;
    vars(1,1:end) = vars_t;
    stds(1,1:end) = stds_t;
    
    for d = 2:files_num
        [header, means(d,1:end), vars(d,1:end), stds(d,1:end)] = single_csv_data(files{d}, num(d,3), num(d,4),properties{1}, 50000, columns);
    end
     means(1:end,1)= num(1:end,1);
     vars(1:end,1)= num(1:end,1);
     stds(1:end,1)= num(1:end,1);

     col = cols;
        

        %axis([0,0,tand(16),10]) ;
        
   %%%  flux plots   %%%%
%         subplot(2,1,1);
% hold on;        
%         %[roer_f, roer_goodness] = fit(means(1:end,1), means(1:end, col),'(a*tand(x))/(1-(tand(x)/b)^2)','lower', [0,0.26]', 'upper', [500,0.4]);
%         [roer_f, roer_goodness] = fit(means(1:end,1), means(1:end, col),'(a*tand(x))/(1-(b*tand(x))*cos(tan(x))^2)');
%        % [exp_f, exp_goodness] = fit(means(1:end,1), means(1:end, col),'exp1');
%         roer_plot = plot(roer_f, means(1:end,1), means(1:end, col));
%        % exp_plot = plot(exp_f, means(1:end,1), means(1:end, col));
%         
%         set(roer_plot(2), 'color', 'r');
        if r_knick > 2
           % [roer1_f, roer1_goodness] = fit(means(1:r_knick,1), means(1:r_knick, col),'(a*tand(x))/(1-(tand(x)/b)^2)', 'lower', [0, 0.26], 'upper', [500, 0.3]);
           [roer1_f, roer1_goodness] = fit(means(1:r_knick,1), means(1:r_knick, col),'(a*tand(x))/(1-(b*tand(x))*cos(tan(x))^2)');
       %      roer1_plot = plot(roer1_f, means(1:r_knick,1), means(1:r_knick, col));

          %  [exp1_f, exp1_goodness] = fit(means(1:r_knick,1), means(1:r_knick, col),'exp1');
            %exp1_plot = plot(exp1_f, means(1:r_knick,1), means(1:r_knick, col));
            
          %  set(roer1_plot(2), 'color', 'k');
        else
           % [roer1_f, roer1_goodness] = fit((1:10).', (1:10).','(a*tand(x))/(1-(tand(x)/b)^2)', 'lower', [0, 0.26], 'upper', [500, 0.3]);
           [roer1_f, roer1_goodness] = fit(means(1:r_knick,1), means(1:r_knick, col),'(a*tand(x))/(1-(b*tand(x))*cos(tan(x))^2)');
            %,'startpoint', [1e-7,means(knick,1)*0.9999], 'lower', [0,means(knick,1)*0.9], 'upper', [1, means(knick,1)*1.1]);
          
           
        end
        
        if s_knick < numel(files)
           
            standa = strcat(num2str(2 * rho * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_d) * ', num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *b)');
            [st2_f, st2_goodness] = fit(means(s_knick:end,1), means(s_knick:end, col),standa);
            %st2_plot = plot(st2_f, means(s_knick:end,1), means(s_knick:end, col));
            %set(st2_plot(2), 'color', 'b');
        else
            standa = strcat(num2str(2 * rho * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_d) * ', num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *b)');
            [st2_f, st2_goodness] = fit((1:10).',(1:10).',standa);
        end
%         
%         
%         
%      e = errorbar(means(1:end,1),means(1:end, col),stds(1:end, col), '.'); 
%      e.LineStyle = 'none';
%         title('Total Flux');
%         xlabel('slope angle');  
%         ylabel('average total flux');
%         legend('location', 'northwest');
%         
%         text((means(1,1)+means(floor(files_num/2),1))/4,means(end, col),...
%           {strcat('Low Roering (k): K=',num2str(roer1_f.a), ' Sc=', num2str(atand(roer1_f.b)),' rsquare=', num2str(roer1_goodness.rsquare)),...                 
%                     strcat('Standa (b): m_d=',num2str(st2_f.m_d), ' beta=', num2str(st2_f.b),' rsquare=', num2str(st2_goodness.rsquare)),...
%                     strcat('Roering (r): K=',num2str(roer_f.a), ' Sc=', num2str(atand(roer_f.b)),' rsquare=', num2str(roer_goodness.rsquare))},...
%           'HorizontalAlignment','left',...
%           'VerticalAlignment','top',...
%           'BackgroundColor',[1 1 1],...
%           'FontSize',10);
%      
%     legend('off');
%       xlswrite(output_file, header, strcat('means_', properties{1}), 'A1');
%     xlswrite(output_file, header, strcat('stds_', properties{1}), 'A1');
%     xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{1}), 'A2');
%     xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{1}), 'A2');
%       

%%%%%%%%%   end flux   %%%%%

%%%     start velocity   %%%
       [header, means_t, vars_t, stds_t] = single_csv_data(files{1}, num(1,3), num(1,4),properties{2},1, columns);
    cols = numel(header);
    means= zeros(files_num, cols);
    vars = zeros(files_num, cols);
    stds = zeros(files_num, cols);
    means(1,1:end) = means_t;
    vars(1,1:end) = vars_t;
    stds(1,1:end) = stds_t;
    for d = 2:files_num
        [header, means(d,1:end), vars(d,1:end), stds(d,1:end)] = single_csv_data(files{d}, num(d,3), num(d,4),properties{2},1,columns);
    end
     means(1:end,1)= num(1:end,1);
     vars(1:end,1)= num(1:end,1);
     stds(1:end,1)= num(1:end,1);

    
      
      hold on;
      
          
   % with flux subplot(2,1,2);
  
    
      xlabel(strcat('mean layer velocity'));
      ylabel('depth layer');
         
    
       furb_vel = 'u*log((0.5-x)/10)';
% before transition angle
      for theta=1:s_knick
          i_fit = fit ((1:cols - 4).',means(theta, 3:cols - 2).', 'exp1');
         % [i_fit,  i_good] = fit ((1:cols - 4).',means(theta, 3:cols - 2).', furb_vel);
          i_plot = plot(i_fit, (1:cols - 4).',means(theta, 3:cols - 2).');
          set(i_plot(2), 'color', colors(theta,1:3));
      end
      
% after transition angle
depth = (h/(cols - 4): h/(cols - 4): h);

standa_v_mb = strcat(num2str(2*rho*g),'* (',num2str(sind(theta)),' - m_d*', num2str(cosd(theta)), ...
                            ')/(b* ', num2str(3*d*sqrt(rho_p/rho*g*cosd(theta))),') * (',num2str(h^1.5),' - (', num2str(h),'-x).^1.5) ');
                        
                       
                        
      for theta=s_knick:files_num-1
         % b0 = sqrt(h)*rho*g*(sind(theta) - num2str(st2_f.m_d)*cosd(theta))/ ... 
            %            (num2str(st2_f.b)*d*sqrt(rho_p*rho*g*cosd(theta)));
   
            
          %i_fit = fit ((1:cols - 4).',means(theta, 3:cols - 2).', standa_v_mb);
       %   i_plot = plot(i_fit, (1:cols - 4).',means(theta, 3:cols - 2).');
       
        %  standa_full = (2*rho*g)* (sind(theta) - st2_f.m_d*cosd(theta)) ...
           %                ./(st2_f.b* 3*d*sqrt(rho_p/rho*g*cosd(theta))) * (h^1.5) - ( h-x).^1.5;
           %  [i_fit, i_good] = fit ((1:cols - 4).',means(theta, 3:cols - 2).', standa_v_mb);
           
           % TODO - read params from output file
          i_plot=plot((1:cols - 4).',  (2*rho*g)* (sind(theta) - st2_f.m_d*cosd(theta)) ...
                         ./(st2_f.b* 3*d*sqrt(rho_p/rho*g*cosd(theta))) * (h^1.5) - ( h-depth.').^1.5);
         %i_plot=plot((1:cols - 4).', i_fit)
         
          %set(i_plot(2), 'color', colors(theta,1:3));
      end
%  text((means(1,1)+means(floor(files_num/2),1))/4,means(end, col),...
%           {strcat('exp  a=',num2str(roer1_f.a), ' Sc=', num2str(atand(roer1_f.b)),' rsquare=', num2str(roer1_goodness.rsquare)),...                 
%                     strcat('Standa (b): m_d=',num2str(st2_f.m_d), ' beta=', num2str(st2_f.b),' rsquare=', num2str(st2_goodness.rsquare)),...
%                     strcat('Roering (r): K=',num2str(roer_f.a), ' Sc=', num2str(atand(roer_f.b)),' rsquare=', num2str(roer_goodness.rsquare))},...
%           'HorizontalAlignment','left',...
%           'VerticalAlignment','top',...
%           'BackgroundColor',[1 1 1],...
%           'FontSize',10);
      
   lh = legend(arrayfun(@num2str, num(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);     
        saveas(gcf, strcat(strrep(output_fig, '.','-'), '_basic'), 'jpg');

   
