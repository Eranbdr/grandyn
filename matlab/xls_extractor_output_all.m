function [] = xls_extractor_output_all()
warning('off','all');

       set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
    colors = get(0, 'DefaultAxesColorOrder');
    vel_property = 'velocity_profile';
    flux_property = 'flux';
    angles = {'_4', '_7', '_10', '_12', '_13', '_13.5', '_14', '_14.5' '_15', '_15.5', '_16'};
    %angles = {'_13'};%, '_15.5', '_16'};
    wls ={{5,  8, 10, 12, 15, 18, 20} { 10, 12, 15}, {10, 12, 15},{ 10, 12}};
    %knicks ={{'18','18','18','18','18','18', '18'}, { '18','18','18'}, {'18','18', '16_15'},{ '16_15', '16_15'}};
    knicks ={{9,9,9,9,9,9,9}, { 9,9,9}, {9,9,8},{ 8,8}};
    %wls ={5, 8};%, 10};
    amps = {11, 19, 25, 31};
     
    %angs=[4,4,7,7,10,10,12,12,13,13,13.5,13.5,14,14,14.5,14.5,15,15,15.5,15.5,16,16];
    %angs_step = 2;
    angs=[4,7,10,12,13,13.5,14,14.5,15,15.5,16];
    angs_step = 1;
    data = zeros(numel(angs),12);
    stds = zeros(numel(angs),12);
    single_profile = zeros(numel(angs),10);
    single_profile_mean = zeros(numel(angs)/2,10);
    sys ='_S2';
    out_col = 'a'; % start at b
    index = 1;
    output_file = strcat('output_all',sys, '.xlsx');
    sheet_name = sys;
    
    % A-f name 	mag dur "A^2/f"	trans angle	ave vel transition	total flux transition
    xlswrite(output_file, {'Name','Magnitude', 'Duration', 'A^2/f', 'trans angle', 'velocity(trans)', 'flux(trans)'},sheet_name, 'a1');
    % H-M  K	Sc		m_d	beta	R rsquare		ST rsquare		angles: total flux
    %xlswrite(output_file, {'K','Sc', 'm_d', 'beta', 'Roer R^2', 'Standa R^2'},sheet_name, 'h1');
    % N-O <angles: total flux>, <angle;a_vel;b_vel;rsquare>
    xlswrite(output_file, {'<angles:total flux>', '<angle, a, b>'},sheet_name, 'n1');
    
    for amp = 1:numel(amps)
        for wl = 1:numel(wls{amp})
            
            % construct name and read file
            file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'gy_', num2str(amps{amp}), 'gx_', num2str(wls{amp}{wl}));
            %sheet_name = strcat('random', sys,  '_', num2str(amps{amp}), 'g_', num2str(wls{amp}{wl}), '_WL');
            [fits, h_fits] = xlsread(strcat(file_name, '000WL_1_knick_',num2str(knicks{amp}{wl}),'_output.xlsx'), 'fits');
            [total_flux_data(1:numel(angs) ,index)] = xlsread(strcat(file_name, '000WL_1_knick_',num2str(knicks{amp}{wl}),'_output.xlsx'), strcat('means_', flux_property), 'm:m');
            [ave_vel_data(1:numel(angs) ,index)] = xlsread(strcat(file_name, '000WL_1_knick_',num2str(knicks{amp}{wl}),'_output.xlsx'), strcat('means_', vel_property), 'm:m');
            %[stds(1:numel(angs) ,index)] = xlsread(strcat(file_name, '000WL_1_knick_',num2str(knicks{amp}{wl}),'_output.xlsx'), strcat('stds_', property), 'm:m');
            [single_profile(1:numel(angs) ,1:10)] = xlsread(strcat(file_name, '000WL_1_knick_',num2str(knicks{amp}{wl}),'_output.xlsx'), strcat('means_', property), 'c:l');
            trans_angle_index = str2num(num2str(knicks{amp}{wl}));
         %   figure('name', file_name);
            %axes1 = subplot(2,1,1);
            
         %   view(axes1,[90 -90]);
          %  hold(axes1,'all');
          
            % single profile (output file)
            for idx = 1:numel(angs) 
                
                %per angle output
                [idx_fit, idx_goodness] = fit ((0:-1:-9).',single_profile_mean(idx, 1:10).', 'exp1');
                
                 xlswrite(output_file, {strcat('<',angs(idx),';',num2str(idx_fit.a),';', num2str(idx_fit.b),';', num2str(idx_goodness.rsquare),'>')} , sheet_name, strcat('n', num2str(1+idx)));
                 xlswrite(output_file, {strcat('<',angs(idx),';',num2str(total_flux_data(idx, index)),'>')} , sheet_name, strcat('o', num2str(1+idx)));
            end
            
            % write to ouput file
              xlswrite(output_file,  {amps{amp}, wls{amp}{wl}, amps{amp}^2* wls{amp}{wl}, angs(trans_angle_index)}, sheet_name, strcat('b', num2str(index + 1)));
              xlswrite(output_file,  { ave_vel_data(trans_angle_index), total_flux_data(trans_angle_index)}, sheet_name, strcat('f', num2str(index + 1)));
              xlswrite(output_file,  fits, sheet_name, strcat('h', num2str(index + 1)));


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

