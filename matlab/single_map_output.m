function [] = single_map_output()
    warning('off','all');
     sys ='_L1';
    wl =12000;
    %wls ={{} { }, {10000, 12000, 15000},{ 10000, 12000}};
    %knicks ={{'18','18','18','18','18','18', '18'}, { '18','18','18'}, {'18','18', '16_15'},{ '16_15', '16_15'}};
     knick =9;
    %wls ={5, 8};%, 10};
    amp = 19;
    output_line = 6;
     
    %angs=[4,4,7,7,10,10,12,12,13,13,13.5,13.5,14,14,14.5,14.5,15,15,15.5,15.5,16,16];
    %angs_step = 2;
    
   
    coeffs = strcat('coeffs', sys);
    flux = strcat('flux', sys);
    flux_std = strcat('flux_std', sys);
    output_file = '../matlab_2/output/final/final_tables.xlsx';
    xlswrite(output_file, {'Magnitude', 'Duration', 'K', 'std K', 'Sc', 'std Sc', 'Furbish R^2', 'mu_d', 'std mu_d', 'beta', 'std beta', 'Parez R^2', 'transition angle'}, coeffs, 'A1');
        
    %file_name = strcat('random', sys,  '_', num2str(amp), 'gy_', num2str(amp), 'gx_', num2str(wl), 'WL_3_output');
   
    file_name = strcat('random_restart', sys,  '_', num2str(amp), 'g_',  num2str(wl), 'WL_4_output');
   %file_name = strcat('random', sys,  '_', num2str(amp), 'g_',  num2str(wl), 'WL_4_output');

       %[trans_angle, trans_velocity, trans_flux, K, Sc, Rr2, Kf, N, Fr2, mu_d, beta, Sr2, flux_pairs1, flux_pairs2] = xls_extract_gather_output(file_name, knicks{amp}{wl});
%        plot_from_output_curr('output/',file_name, knicks{amp}{wl}, 10, 0);
        [Kf, stdK, Sc, stdSc, Fr2, mu_d, stdM, beta, stdB, Sr2, fluxes, fluxes_stds, angles, trans_angle]= coeffs_from_output_curr('../matlab_2/output/',file_name, knick, 20);

        xlswrite(output_file, {'Magnitude', 'Duration'}, flux, 'A1');
    xlswrite(output_file, angles.', flux, 'C1');
    xlswrite(output_file, {'Magnitude', 'Duration'}, flux_std, 'A1');
    xlswrite(output_file, angles.', flux_std, 'C1');

        xlswrite(output_file, {num2str(amp) num2str(wl) num2str(Kf) num2str(stdK) num2str(Sc) num2str(stdSc) num2str(Fr2) num2str(mu_d) num2str(stdM) num2str(beta) num2str(stdB) num2str(Sr2) trans_angle}, coeffs, strcat('A', num2str(output_line)));
        xlswrite(output_file, {num2str(amp) num2str(wl)}, flux, strcat('A', num2str(output_line)));
        xlswrite(output_file, {num2str(amp) num2str(wl)}, flux_std, strcat('A', num2str(output_line)));
        xlswrite(output_file, fluxes.', flux, strcat('C', num2str(output_line)));
        xlswrite(output_file, fluxes_stds.', flux_std, strcat('C', num2str(output_line)));


      
end