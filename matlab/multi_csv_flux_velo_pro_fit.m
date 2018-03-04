function[] = multi_csv_flux_velo_pro_fit(sub_name)

% read several CSVs and present means and vars 
% input is slope, file, SS_step
warning('off','all');
solid_fraction = 0.8;%965/(34*34);
g = 10;
d = 0.001;
rho_p = 2500;
rho = rho_p * solid_fraction;
h = 34*0.001;

columns = 23;% S 13, M 23

r_knick = 6; % from lowest slope
s_knick = 5;
set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
colors = get(0, 'DefaultAxesColorOrder');

properties = {'flux', 'velocity_profile'};

[num, files] = xlsread(strcat('C:\Users\EranBD\Documents\thesis\matlab\', sub_name, '_input.csv'));
output_file = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', sub_name, '_output.xlsx');
output_fig = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', sub_name, '_knick_', num2str(r_knick));
% [num, files] = xlsread(strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\dups\', sub_name, '_input.csv'));
% output_file = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\dups\', sub_name, '_output.xlsx');
% output_fig = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\dups\', sub_name);
delete output_file;
files_num = numel(files);
figure1 = figure('name', strcat(sub_name, '_knick_', num2str(r_knick)));


    
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
      xlswrite(output_file, header, strcat('means_', properties{1}), 'A1');
    xlswrite(output_file, header, strcat('stds_', properties{1}), 'A1');
    xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{1}), 'A2');
    xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{1}), 'A2');
    %  return;
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

     xlswrite(output_file, header, strcat('means_', properties{2}), 'A1');
     xlswrite(output_file, header, strcat('stds_', properties{2}), 'A1');
     xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{2}), 'A2');
     xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{2}), 'A2');
      

      
          
    subplot(2,1,2);
      hold on;  
    
      xlabel(strcat('mean layer velocity'));
      ylabel('depth layer');
         
    
      
% before transition angle
      for theta=1:s_knick
          i_fit = fit ((1:cols - 4).',means(theta, 3:cols - 2).', 'exp1');
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
       
       %   standa_full = (2*rho*g)* (sind(theta) - st2_f.m_d*cosd(theta)) ...
          %                 /(st2_f.b* 3*d*sqrt(rho_p/rho*g*cosd(theta))) * (h^1.5) - ( h-x).^1.5;
          i_plot=plot((1:cols - 4).',  (2*rho*g)* (sind(theta) - st2_f.m_d*cosd(theta)) ...
                           ./(st2_f.b* 3*d*sqrt(rho_p/rho*g*cosd(theta))) * (h^1.5) - ( h-depth.').^1.5);
          set(i_plot(2), 'color', colors(theta,1:3));
      end
 text((means(1,1)+means(floor(files_num/2),1))/4,means(end, col),...
          {strcat('exp  a=',num2str(roer1_f.a), ' Sc=', num2str(atand(roer1_f.b)),' rsquare=', num2str(roer1_goodness.rsquare)),...                 
                    strcat('Standa (b): m_d=',num2str(st2_f.m_d), ' beta=', num2str(st2_f.b),' rsquare=', num2str(st2_goodness.rsquare)),...
                    strcat('Roering (r): K=',num2str(roer_f.a), ' Sc=', num2str(atand(roer_f.b)),' rsquare=', num2str(roer_goodness.rsquare))},...
          'HorizontalAlignment','left',...
          'VerticalAlignment','top',...
          'BackgroundColor',[1 1 1],...
          'FontSize',10);
      
   lh = legend(arrayfun(@num2str, num(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);     
        saveas(gcf, strcat(strrep(output_fig, '.','-'), '_basic'), 'jpg');

     xlswrite(output_file, header, strcat('means_', properties{2}), 'A1');
     xlswrite(output_file, header, strcat('stds_', properties{2}), 'A1');
     xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{2}), 'A2');
     xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{2}), 'A2');
     xlswrite(output_file, {'' 'K' 'Sc'  'm_d' 'beta' 'a' 'b' 'R rsquare'  'ST rsquare'}, 'fits', 'A1');
     xlswrite(output_file, {'' num2str(roer1_f.a) num2str(roer1_f.b)  num2str(st2_f.m_d) num2str(st2_f.b) num2str(roer1_goodness.rsquare) num2str(st2_goodness.rsquare)}, 'fits', 'A2');
   
