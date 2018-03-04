function [header, means, vars, stds] = single_csv_data(file_name, SS_first_step, SS_last_Step,property, factor, cols_num)
    [csv_data_orig, header] = xlsread(strcat('../matlab/',property, '_', file_name, '.csv'));%,strcat('A:',char('a' + cols_num)));
   
    SS_step_index = find(csv_data_orig() == SS_first_step, 1);
    means = zeros(cols_num,1);
    vars = zeros(cols_num,1);
    stds = zeros(cols_num,1);
    if factor == 1
        dt_f = 1;
    else
        dt_f = (csv_data_orig(3,1) - csv_data_orig(2,1))/factor;
    end
     [rows, cols] = size(csv_data_orig);
     last_index = find(csv_data_orig() == SS_last_Step, 1);
     if  numel(last_index) == 0
          last_index = rows;
     end
   
    csv_data = csv_data_orig(SS_step_index:last_index, 1:cols_num);   
    csv_data(isnan(csv_data)) = 0 ;
     if sum(csv_data(:) == 0) > 1000
         display(strcat('too many NaNs, ', file_name));
         if strcmp(property, 'velocity_profile') == 0
             if sum(csv_data(1:end,cols_num) == 0) > 100
                 display(strcat('even for total flux ', file_name));
                 return;
             else
                 display(strcat('but ok for total flux: ', num2str( sum(csv_data(1:end,cols_num) == 0))));
             end
         else
             display(strcat('in velocity: ', num2str( sum(csv_data(1:end,cols_num) == 0))));
             return;
         end
     end
    for curr = 2:cols_num
        means(curr) = mean(csv_data(1:end , curr)./dt_f);
        vars(curr) = var(csv_data(1:end , curr)./dt_f);
        stds(curr) = std(csv_data(1:end , curr)./dt_f);
    end
    
    % write to ouput file
    %output_file = strcat('single_', file_name, '_', property, '.xlsx');
     %xlswrite(output_file, header, strcat('means_', property), 'A1');
     %xlswrite(output_file, header, strcat('stds_', property), 'A1');
     %xlswrite(output_file, means(1:end,1:end).', strcat('means_', property), 'A2');
     %xlswrite(output_file, stds(1:end,1:end).', strcat('stds_', property), 'A2');
end

