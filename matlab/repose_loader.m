function [] = repose_loader()

fromS = [000105200000 000166100000 000203750000 000429400000 000173650000]/50000+15000;
fromM = [000431150000 000300100000 000737900000 000639550000 000302700000]/50000+15000;
fromL = [000630700000 001113750000 000292800000 000254400000 001128350000]/50000+15000;

from = fromM;
    for i = 1:5
        file_name = strcat('C:\Users\EranBD\Documents\thesis\repose\angle_of_repose_M_', ...
        num2str(i), '_restart\flux_angle_of_repose_M_', num2str(i),'_restart.csv');
        [csv_data, header] = xlsread(file_name);%,strcat('A:',char('a' + cols_num)));

        display(num2str(mean(csv_data(from(i):end , 4))));
    end
    
end

