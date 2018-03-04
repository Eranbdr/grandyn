function[] =multi_csv_totflux_velpro_no_plot(sub_name)
% read several CSVs and present means and vars 
% input is slope, file, SS_step
warning('off','all');
%sub_name = 'random_S2_11gy_11gx_12000WL_1';
columns = 23;% S 13, M 23

properties = {'flux', 'velocity_profile'};

[num, files] = xlsread(strcat('C:\Users\EranBD\Documents\thesis\matlab\', sub_name, '_input.csv'));
output_file = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', sub_name,'_output_2d.xlsx');
delete output_file;
files_num = numel(files);


[header, means_t, vars_t, stds_t] = single_csv_data(files{1}, num(1,3), num(1,4), properties{1}, 50000, columns);
 
means= zeros(files_num, columns);
vars = zeros(files_num, columns);
stds = zeros(files_num, columns);
means(1,1:end) = means_t;
vars(1,1:end) = vars_t;
stds(1,1:end) = stds_t;

for d = 2:files_num
    [header, means(d,1:end), vars(d,1:end), stds(d,1:end)] = single_csv_data(files{d}, num(d,3), num(d,4),properties{1}, 50000, columns);
end
means(1:end,1)= num(1:end,1);
vars(1:end,1)= num(1:end,1);
stds(1:end,1)= num(1:end,1);

xlswrite(output_file, header, strcat('means_', properties{1}), 'A1');
xlswrite(output_file, header, strcat('stds_', properties{1}), 'A1');
xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{1}), 'A2');
xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{1}), 'A2');
   

[header, means_t, vars_t, stds_t] = single_csv_data(files{1}, num(1,3), num(1,4),properties{2},1, columns);

means= zeros(files_num, columns);
vars = zeros(files_num, columns);
stds = zeros(files_num, columns);
means(1,1:end) = means_t;
vars(1,1:end) = vars_t;
stds(1,1:end) = stds_t;
for d = 2:files_num
    [header, means(d,1:end), vars(d,1:end), stds(d,1:end)] = single_csv_data(files{d}, num(d,3), num(d,4),properties{2},1,columns);
end
means(1:end,1)= num(1:end,1);
vars(1:end,1)= num(1:end,1);
stds(1:end ,1)= num(1:end,1);

xlswrite(output_file, header, strcat('means_', properties{2}), 'A1');
xlswrite(output_file, header, strcat('stds_', properties{2}), 'A1');
xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{2}), 'A2');
xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{2}), 'A2');
    
