function [] = xls_extractor_output()
warning('off','all');

       set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
    colors = get(0, 'DefaultAxesColorOrder');
    property = 'velocity_profile';
    %property = 'vol_flux';
    angles = {'_4', '_7', '_10', '_12', '_13', '_13.5', '_14', '_14.5' '_15', '_15.5', '_16'};
    %angles = {'_13'};%, '_15.5', '_16'};
    wls ={{5,  8, 10, 12, 15, 18, 20} { 10, 12, 15}, {10, 12, 15},{ 10, 12}};
    knicks ={{'18','18','18','18','18','18', '18'}, { '18','18','18'}, {'18','18', '16_15'},{ '16_15', '16_15'}};
    %wls ={5, 8};%, 10};
    amps = {11, 19, 25, 31};
    data = zeros(22,12);
    stds = zeros(22,12);
    angs=[4,4,7,7,10,10,12,12,13,13,13.5,13.5,14,14,14.5,14.5,15,15,15.5,15.5,16,16];
    
    single_profile = zeros(numel(angs),10);
    single_profile_mean = zeros(numel(angs)/2,10);
    sys ='_S12';
    out_col = 'a'; % start at b
    index = 1;
    output_file = strcat('outputs',sys, '_', property, '.xlsx');
    A_2WL = zeros(numel(v_sq));
    for amp = 1:numel(amps)
        for wl = 1:numel(wls{amp})

            file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'gy_', num2str(amps{amp}), 'gx_', num2str(wls{amp}{wl}));
            sheet_name = strcat('random', sys,  '_', num2str(amps{amp}), 'g_', num2str(wls{amp}{wl}), '_WL');
            [fits, h_fits] = xlsread(strcat(file_name, '000WL_1_knick_',knicks{amp}{wl},'_output.xlsx'), 'fits');
            [data(1:numel(angs) ,index)] = xlsread(strcat(file_name, '000WL_1_knick_',knicks{amp}{wl},'_output.xlsx'), strcat('means_', property), 'm:m');
            [stds(1:numel(angs) ,index)] = xlsread(strcat(file_name, '000WL_1_knick_',knicks{amp}{wl},'_output.xlsx'), strcat('stds_', property), 'm:m');

            [single_profile(1:numel(angs) ,1:10)] = xlsread(strcat(file_name, '000WL_1_knick_',knicks{amp}{wl},'_output.xlsx'), strcat('means_', property), 'c:l');
         %   figure('name', file_name);
            %axes1 = subplot(2,1,1);
            
         %   view(axes1,[90 -90]);
          %  hold(axes1,'all');
            for idx = 1:2:numel(angs) 
                single_profile_mean(ceil((idx+1)/2), 1:end) = mean(single_profile(idx:idx+1, 1:end));
                [idx_fit, idx_goodness] = fit ((1:10).',single_profile_mean(ceil((idx+1)/2), 1:10).', 'exp1');
               % idx_plot = plot(idx_fit, (1:10).',single_profile_mean(ceil((idx+1)/2), 1:10).');
           %     set(idx_plot(2), 'color', colors(ceil((idx+1)/2),1:3));
                 xlswrite('velocity_profile_by_angle.xlsx', {'run','vel_pro_a', 'vel_pro_b', 'vel_pro_R^2'},num2str(angs(idx)), 'a1');
                 xlswrite('velocity_profile_by_angle.xlsx',  {sheet_name}, num2str(angs(idx)),strcat('A', num2str(1+index)));
                 %xlswrite('velocity_profile.xlsx',num2str(angs(idx)),  {strcat(sheet_name, '_WL')}, strcat('b', num2str(ceil(1+(idx+1)/2))));
                 xlswrite('velocity_profile_by_angle.xlsx', {num2str(idx_fit.a), num2str(idx_fit.b), num2str(idx_goodness.rsquare)},  num2str(angs(idx)), strcat('b',num2str(1+index)));
            end
            
            
           % ylabel('layer');
%             legend('off');
%             axes2 = subplot(2,1,2);
%             
%             view(axes2,[90 -90]);
%             hold(axes2,'all');
%             for idx = 1:2:numel(angs) 
%                       idx1_plot = plot((1:10).',(single_profile_mean(ceil((idx+1)/2), 1:10)./single_profile_mean(ceil((idx+1)/2),10)).');
%                       set(idx1_plot, 'color', colors(ceil((idx+1)/2),1:3));
%             end
%             ylabel('normalized velocity');
%             xlabel('layer');
%             lh = legend(angles, 'location', 'southeast', 'FontSize',7);
%             saveas(gcf, strcat('vel_pro_',file_name));
            % write to ouput file
%             xlswrite(output_file, {'amp', 'WL', 'WLA^2','v_sq'}, 'fits', 'A1');
%              xlswrite(output_file, h_fits, 'fits', 'E1');
%              xlswrite(output_file,  {amps{amp}, wls{amp}{wl}, amps{amp}^2* wls{amp}{wl}, v_sq(index)}, 'fits', strcat('A', num2str(index + 1)));
%              A_2WL(index) = amps{amp}^2* wls{amp}{wl};
%              xlswrite(output_file, fits, 'fits', strcat('E', num2str(index + 1)));


             index = index + 1;

        end
    end
      
    means= zeros(ceil(numel(angs)/2), index - 1);
    xlswrite(output_file, v_sq, property, 'B1');
    xlswrite(output_file, {'P1','P2', 'Rsq'}, 'fits_b', 'B1');
    fits = cell(1,numel(angs)/2);
    goodness = cell(1,numel(angs)/2);
    figure;
    hold on;
    for idx = 1:2:numel(angs) 
        means(ceil((idx+1)/2), 1:end) = mean(data(idx:idx+1, 1:end));
        xlswrite(output_file, angs(idx), property, strcat('A', num2str(ceil((idx+1)/2) +1)));
        xlswrite(output_file, means(ceil((idx+1)/2), 1:end),property, strcat('B', num2str(ceil((idx+1)/2) + 1)));
       [fits{ceil((idx+1)/2)},goodness{ceil((idx+1)/2)}] = fit(v_sq(1:index-1).',means(ceil((idx+1)/2), 1:end).','poly1');
       p = plot(fits{ceil((idx+1)/2)},v_sq(1:index-1),means(ceil((idx+1)/2),1:end));
       xlswrite(output_file, angs(idx), 'fits_b', strcat('A', num2str(ceil((idx+1)/2) +1)));
        xlswrite(output_file, {num2str(fits{ceil((idx+1)/2)}.p1), num2str(fits{ceil((idx+1)/2)}.p2), num2str(goodness{ceil((idx+1)/2)}.rsquare)}, 'fits_b', strcat('B', num2str(ceil((idx+1)/2) + 1)));
       set(p(2), 'color', colors(ceil((idx+1)/2),1:3));
       e = errorbar(v_sq(1:index-1),means(ceil((idx+1)/2), 1:end), std(data(idx:idx+1, 1:end)), '.', 'color', colors(ceil((idx+1)/2),1:3));
       e.LineStyle = 'none';
       e.markeredgecolor=colors(ceil((idx+1)/2));
       e.markerfacecolor=colors(ceil((idx+1)/2));
    end
    
    xlabel('V^2');
    ylabel(property);
 legend(arrayfun(@num2str, angs(end:-2:1,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
 %   legend(arrayfun(@num2str, angs(1:2:end), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
end

