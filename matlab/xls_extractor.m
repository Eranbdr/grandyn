function [header, means, vars, stds] = xls_extractor()

    %property = 'velocity_profile';
    property = 'vol_flux';
    sheet = 'means_vol_flux';
    
    
    
    angles = {'_4', '_7', '_10', '_12', '_13', '_13.5', '_14', '_14.5' , '_15', '_15.5', '_16'};
    %angles = {'_13', '_13.5', '_14.5'};
    wl ={'_10', '_12', '_15', '_18', '_20'};
    wl ={'_10', '_12', '_15'};
    amps = {'_11g', '_19g', '_25g', '_31g'};
    sys ='_S1';
    out_col = 'a'; % start at b
    first_index = 10000;
    last_index = 40000;
    first_index = 2;
    last_index = 9;
    output_file = strcat('AMP_alls-angles',sys, '.xlsx');
    for amp = 1:numel(amps)
        for i = 1:numel(wl)
            for angle = 1:numel(angles)
                file_name = strcat('random', sys, amps{amp}, 'y', amps{amp}, 'x', wl{i});
                [csv_data, header] = xlsread(strcat(file_name, '000WL_1_knick_9_output.xlsx'), sheet,  'm:m');

                [rows, cols] = size(csv_data);
                
              %  last_index = min(last_index, rows);
%                 for curr = 1:1
%                     means(curr) = mean(csv_data(first_index:last_index , curr));
%                     vars(curr) = var(csv_data(first_index:last_index , curr));
%                     stds(curr) = std(csv_data(first_index:last_index , curr));
%                 end

                % write to ouput file
                 %xlswrite(output_file, {strcat(sys,angles{a},amp)}, strcat(char(a+out_col),'1'));
                 xlswrite(output_file, {strcat(sys,angles{angle})}, property, strcat(char(angle+out_col),'1'));
                 %xlswrite(output_file, csv_data.', 'g',  strcat(char(a + out_col), '2'));
                 xlswrite(output_file, csv_data.', property,  strcat('b', num2str((amp-1)*numel(wl)+i+1)));
            end
         
        end
        xlswrite(output_file, strcat(amps{amp}, wl).', property, strcat('A', num2str((amp-1)*numel(wl)+2)));
    end
end

