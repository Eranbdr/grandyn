function [] = generic_map_output()
    warning('off','all');
     sys ='_S1';
    wls ={{ 5000, 8000, 10000, 12000, 15000, 18000, 20000} { 10000, 12000, 15000}, {10000, 12000, 15000},{ 10000, 12000}};
    %wls ={{} { }, {10000, 12000, 15000},{ 10000, 12000}};
    %knicks ={{'18','18','18','18','18','18', '18'}, { '18','18','18'}, {'18','18', '16_15'},{ '16_15', '16_15'}};
    
    if strcmp(sys, '_S1')
        %S1 11 10 missing angle 14.5
         knicks ={{9,9,8,9,9,9,9}, { 9,9,9}, {9,9,9},{ 9,8}};
    elseif strcmp(sys, '_S2')
            % no angle 13 except S2 11 15
            knicks ={{8,8,8,8,9,8,8}, { 8,8,8}, {8,8,8},{ 8,7}};
    elseif strcmp(sys, '_S3')
            
            knicks ={{9,9,9,9,9,9,9}, { 9,9,9}, {9,9,9},{ 8,9}};
    end
    %wls ={5, 8};%, 10};
    amps = {11, 19, 25, 31};
     
    %angs=[4,4,7,7,10,10,12,12,13,13,13.5,13.5,14,14,14.5,14.5,15,15,15.5,15.5,16,16];
    %angs_step = 2;
    
   
    coeffs = strcat('coeffs', sys);
    flux = strcat('flux', sys);
    flux_std = strcat('flux_std', sys);
    Ek = strcat('Ek', sys);
    ratio = strcat('ratio', sys);
    %output_file = '../matlab_2/output/final/final_tables_Ek2.xlsx';
    output_file = '../matlab_2/ooutput/final_tables_Ek2.xlsx';
    xlswrite(output_file, {'Magnitude', 'Duration', 'K', 'std K', 'Sc', 'std Sc', 'Furbish R^2', 'mu_d', 'std mu_d', 'beta', 'std beta', 'Parez R^2', 'transition angle'}, coeffs, 'A1');
    index = 1;
    for amp = 1:numel(amps)
        for wl = 1:numel(wls{amp})
            index = index + 1;
            % construct name and read file
            
            % gy gx , knick format
            if strcmp(sys, '_S1')
                file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'gy_', num2str(amps{amp}), 'gx_', num2str(wls{amp}{wl}), 'WL_1_output_2d');
                %file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'gy_', num2str(amps{amp}), 'gx_', num2str(wls{amp}{wl}), 'WL_1_knick_', num2str(knicks{amp}{wl}),'_output');
            
            % gy gx , no knick format
            elseif strcmp(sys, '_S2')
                file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'gy_', num2str(amps{amp}), 'gx_', num2str(wls{amp}{wl}), 'WL_1_output_2d');
            
            % g , 2d format
            elseif strcmp(sys, '_S3')
                file_name = strcat('random', sys,  '_', num2str(amps{amp}), 'g_', num2str(wls{amp}{wl}), 'WL_5_output_2d');
            end
           %[trans_angle, trans_velocity, trans_flux, K, Sc, Rr2, Kf, N, Fr2, mu_d, beta, Sr2, flux_pairs1, flux_pairs2] = xls_extract_gather_output(file_name, knicks{amp}{wl});
%        plot_from_output_curr('output/',file_name, knicks{amp}{wl}, 10, 0);
            [Kf, stdK, Sc, stdSc, Fr2, mu_d, stdM, beta, stdB, Sr2, fluxes, fluxes_stds, angles, trans_angle, Eks]= coeffs_from_output_curr('../matlab_2/ooutput/',file_name, knicks{amp}{wl}, 10);
            
            if index == 2
                xlswrite(output_file, {'Magnitude', 'Duration'}, flux, 'A1');
                xlswrite(output_file, angles.', flux, 'C1');
                xlswrite(output_file, {'Magnitude', 'Duration'}, flux_std, 'A1');
                xlswrite(output_file, angles.', flux_std, 'C1');
                xlswrite(output_file, {'Magnitude', 'Duration'}, Ek, 'A1');
                xlswrite(output_file, angles.', Ek, 'C1');
                xlswrite(output_file, {'Magnitude', 'Duration'}, ratio, 'A1');
                xlswrite(output_file, angles.', ratio, 'C1');
            end
            xlswrite(output_file, {num2str(amps{amp}) num2str(wls{amp}{wl}) num2str(Kf) num2str(stdK) num2str(Sc) num2str(stdSc) num2str(Fr2) num2str(mu_d) num2str(stdM) num2str(beta) num2str(stdB) num2str(Sr2) trans_angle}, coeffs, strcat('A', num2str(index)));
            xlswrite(output_file, {num2str(amps{amp}) num2str(wls{amp}{wl})}, flux, strcat('A', num2str(index)));
            xlswrite(output_file, {num2str(amps{amp}) num2str(wls{amp}{wl})}, flux_std, strcat('A', num2str(index)));
            xlswrite(output_file, {num2str(amps{amp}) num2str(wls{amp}{wl})}, Ek, strcat('A', num2str(index)));
            xlswrite(output_file, {num2str(amps{amp}) num2str(wls{amp}{wl})}, ratio, strcat('A', num2str(index)));
            xlswrite(output_file, fluxes.', flux, strcat('C', num2str(index)));
            xlswrite(output_file, fluxes_stds.', flux_std, strcat('C', num2str(index)));
            xlswrite(output_file, Eks.', Ek, strcat('C', num2str(index)));
            xlswrite(output_file, (Eks./fluxes).', ratio, strcat('C', num2str(index)));
            
            
        end
    end
end