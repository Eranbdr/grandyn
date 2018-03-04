function [] = xls_extractor_map_output()
    warning('off','all');

    sheet_name = 'table';
    
    wls ={{5000,  8000, 10000, 12000, 15000, 18000, 20000} { 10000, 12000, 15000}, {10000, 12000, 15000},{ 10000, 12000}};
    %knicks ={{'18','18','18','18','18','18', '18'}, { '18','18','18'}, {'18','18', '16_15'},{ '16_15', '16_15'}};
    knicks ={{9,9,9,9,9,9,9}, { 9,9,9}, {9,9,8},{ 8,8}};
    %wls ={5, 8};%, 10};
    amps = {11, 19, 25, 31};
     
    %angs=[4,4,7,7,10,10,12,12,13,13,13.5,13.5,14,14,14.5,14.5,15,15,15.5,15.5,16,16];
    %angs_step = 2;
    
    sys ='_S12';
    
    
    output_file = strcat('output/output_all',sys, '.xlsx');
    
      
    xlswrite(output_file, {'system'}, sheet_name, 'b1');
    xlswrite(output_file, {'magnitude'}, sheet_name, 'c1');
    xlswrite(output_file, {'duration'}, sheet_name, 'd1');
    xlswrite(output_file, {'transtion_angle'}, sheet_name, 'e1');
    xlswrite(output_file, {'trans_velocity'}, sheet_name, 'f1');
    xlswrite(output_file, {'trans_flux'}, sheet_name, 'g1');
    xlswrite(output_file, {'K'}, sheet_name, 'h1');
    xlswrite(output_file, {'Sc'}, sheet_name, 'i1');
    xlswrite(output_file,{ 'Rr2'}, sheet_name, 'j1');
    xlswrite(output_file, {'Kf'}, sheet_name, 'k1');
    xlswrite(output_file, {'N'}, sheet_name, 'l1');
    xlswrite(output_file, {'Fr2'}, sheet_name, 'm1');
    xlswrite(output_file, {'mu_d'}, sheet_name, 'n1');
    xlswrite(output_file, {'beta'}, sheet_name, 'o1');
    xlswrite(output_file, {'Sr2'}, sheet_name, 'p1');
    xlswrite(output_file, {'flux_pairs'}, sheet_name, 'q1');
    
    index = 2;
    
    for amp = 1:numel(amps)
        for wl = 1:numel(wls{amp})
            
            % construct name and read file
            file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'gy_', num2str(amps{amp}), 'gx_', num2str(wls{amp}{wl}), 'WL_1_knick_', num2str(2*knicks{amp}{wl}),'_output.xlsx');
           [trans_angle, trans_velocity, trans_flux, K, Sc, Rr2, Kf, N, Fr2, mu_d, beta, Sr2, flux_pairs1, flux_pairs2] = xls_extract_gather_output(file_name, knicks{amp}{wl});
            
           xlswrite(output_file, {'S1'}, sheet_name, strcat('b', num2str(index)));
           xlswrite(output_file, {'S2'}, sheet_name, strcat('b', num2str(index+1)));
           xlswrite(output_file, amps{amp}, sheet_name, strcat('c', num2str(index)));
           xlswrite(output_file, amps{amp}, sheet_name, strcat('c', num2str(index+1)));
           xlswrite(output_file, wls{amp}{wl}, sheet_name, strcat('d', num2str(index)));
           xlswrite(output_file, wls{amp}{wl}, sheet_name, strcat('d', num2str(index+1)));
           xlswrite(output_file, {flux_pairs1}, sheet_name, strcat('q', num2str(index)));
           xlswrite(output_file, {flux_pairs2}, sheet_name, strcat('q', num2str(index+1)));
           
           for i = 1:2
                
                
                 xlswrite(output_file, {trans_angle(i)}, sheet_name, strcat('e', num2str(index)));
                 xlswrite(output_file, {trans_velocity(i)}, sheet_name, strcat('f', num2str(index)));
                 xlswrite(output_file, {trans_flux(i)}, sheet_name, strcat('g', num2str(index)));
                 xlswrite(output_file, {K(i)}, sheet_name, strcat('h', num2str(index)));
                 xlswrite(output_file, {Sc(i)}, sheet_name, strcat('i', num2str(index)));
                 xlswrite(output_file, {Rr2(i)}, sheet_name, strcat('j', num2str(index)));
                 xlswrite(output_file, {Kf(i)}, sheet_name, strcat('k', num2str(index)));
                 xlswrite(output_file, {N(i)}, sheet_name, strcat('l', num2str(index)));
                 xlswrite(output_file, {Fr2(i)}, sheet_name, strcat('m', num2str(index)));
                 xlswrite(output_file, {mu_d(i)}, sheet_name, strcat('n', num2str(index)));
                 xlswrite(output_file, {beta(i)}, sheet_name, strcat('o', num2str(index)));
                 xlswrite(output_file, {Sr2(i)}, sheet_name, strcat('p', num2str(index)));
                 
                 index = index + 1;
           end
            

            
        end
    end
end