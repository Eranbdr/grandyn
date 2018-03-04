function[] = multi_csv_velo_pro_fit(sub_name, )



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
delete output_file;
files_num = numel(files);
figure('name', strcat(sub_name, '_knick_', num2str(r_knick)));


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


      %p = plot( fliplr(-(cols-4):-1).',(means(1:end,3:cols - 2)./repmat(means(1:end,cols-2),1,cols - 4)).', '.');
      %p = plot( (1:cols-3).',(means(1:end,3:cols - 1)).', '.');
      %rotate(p, [0 0 1], 270, [0 0 0]);
      %axis([1,10,0,5]);
      %hold on;
      %for i = 1:files_num
          %plot((means(1:2:end, 3:cols - 1)).', 1:cols - 3, '.');
          %plot((means(2:2:end, 3:cols - 1)).', 1:cols - 3, 's','MarkerSize',3);
          
        %  plot((means(1:end, 3:cols - 1)).', 1:cols - 3, '.');
          
          %e = errorbar((1:cols - 3).',(means(i,3:cols - 1)).',stds(i, 3:cols-1).', '.');
          %e.LineStyle = 'none';
          %set(e, 'MarkerFaceColor',colors(i, 1:3));
      %end
      title('Velocity Profile');
      xlabel(strcat('mean layer velocity'));
      ylabel('depth layer');
        lh = legend(arrayfun(@num2str, num(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
      %lh = legend(arrayfun(@num2str, num(1:2:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
     
      hold off;     
     % saveas(gcf, strrep(strcat('postr_',sub_name, '_knick_', num2str(r_knick)), '.','-'), 'jpg');
      %saveas(gcf, strcat('postr_',sub_name, '_knick_', num2str(r_knick)));
     figure1=  figure('name', strcat(sub_name, '_vel_fit'));
      axes1 = axes('Parent',figure1);
view(axes1,[90 -90]);
%hold(axes1,'all');
      for i=1:files_num - 1
          i_fit = fit ((1:cols - 4).',means(i, 3:cols - 2).', 'exp1');
          i_plot = plot(i_fit, (1:cols - 4).',means(i, 3:cols - 2).');
          set(i_plot(2), 'color', colors(i,1:3));
      end
 xlabel('velocity');
 ylabel('layer');
   lh = legend(arrayfun(@num2str, num(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);     
   
   
   
   %%%  normalized figure   %%%%%
%         figure2=  figure('name', strcat(sub_name, '_vel_norm'));
%         title(strcat('raw data: ',strrep(sub_name,'_', ' ')));
%       axes2 = axes('Parent',figure2);
% view(axes2,[90 -90]);
% hold(axes2,'all');
% for i=1:files_num 
%           plot((1:cols - 3).',(means(i, 3:cols - 1)./means(i,cols-1)).');
%           %set(i_plot(2), 'color', colors(i,1:3));
% end
%  ylabel('normalized velocity');
%  xlabel('layer');
%    lh = legend(arrayfun(@num2str, num(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
    %%%%%  normalized    %%%%
    %end

%      xlswrite(output_file, header, strcat('means_', properties{2}), 'A1');
%      xlswrite(output_file, header, strcat('stds_', properties{2}), 'A1');
%      xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{2}), 'A2');
%      xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{2}), 'A2');
%      xlswrite(output_file, {'' 'K' 'Sc' 'I' 'm_2' 'm_s' 'm_d' 'beta' 'R rsquare'  'F2 rsquare' 'ST rsquare'}, 'fits', 'A1');
%     xlswrite(output_file, {'' num2str(roer1_f.a) num2str(roer1_f.b)  num2str(fern2_f.I) num2str(fern2_f.m_2) num2str(fern2_f.m_s) num2str(st2_f.m_d) num2str(st2_f.b) num2str(roer1_goodness.rsquare) num2str(fern2_goodness.rsquare) num2str(st2_goodness.rsquare)}, 'fits', 'A2');
   
