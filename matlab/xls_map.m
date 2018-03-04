function [] = xls_map()

    set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
    colors = get(0, 'DefaultAxesColorOrder');
    property = 'velocity_profile';
    %property = 'vol_flux';
    angles = {'_4', '_7', '_10', '_12', '_13', '_13.5', '_14', '_14.5' '_15', '_15.5', '_16'};
    angles = {'_13'};%, '_15.5', '_16'};
    wls ={{'_5',  '_12', '_18', '_20'}, { '_10', '_12', '_15'}, {'_10', '_12', '_15'},{ '_10', '_12'}};
    %wls ={'_5', '_8'};%, '_10'};
    amps = {'_11g', '_19g', '_25', '_31'};
    sys ='_S12';
    first_index = 90000;
    last_index = 150000;
    first_col = 'm';
    last_col = 'm';
    dt = 10000;
    sheet = '';
    figure('name', strcat(property, angles{1}));
    hold on;
    for amp = 1:numel(amps)
        for angle = 1:numel(angles)
            for wl = 1:numel(wls{amp})
                file_name = strcat(property, '_random', sys, angles{angle}, amps{amp}, 'y', amps{amp}, 'x', wls{amp}{wl}, '000WL_1.csv');
                display (file_name);
               %[p_fft, p_raw  ] = xls_fft_single_fig(file_name, sheet, first_col, last_col, first_index, last_index, dt); 
                %try
                   % set(p_fft(2), 'color', colors(wl + numel(wls{amp-1}), 1:3));
                    %set(p_raw(2), 'color', colors(wl + numel(wls{amp-1}), 1:3));
                %catch ex
                   % set(p_fft(2), 'color', colors(wl), 1:3);
                    %set(p_raw(2), 'color', colors(wl), 1:3);    
                %end
                 xls_fix_nan(file_name, sheet, first_col, last_col, first_index, last_index, dt); 
            end
        end
    end
end

