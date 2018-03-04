function[K, stdK, Sc, stdSc, r_furb, Mu, stdMu,  Beta, stdBeta, r_standa, fluxes, fluxes_stds, angles, trans_angle, Ek] = coeffs_from_output_curr(dir, file_name, transition_angle, num_of_layers, Ek)
%       
warning('off','all');
solid_fraction_S = 0.765;
solid_fraction_ML = 0.78;
solid_fraction = solid_fraction_ML;
g = 2e-10;
d = 1;
total_col = num_of_layers + 3;%to_col + 1;
h = 68.8;
%h = 34.5;
dt1 = 0.07158707624543059;
dt2 = 7.1570724900297727E-002;
dt3 = 0.071561356128674;
dtL1 = 7.16E-02;
dtM1 = 0.0715685656081002;
dt = dtM1;

% time steps and average volume

% skip over step and first column
from_col = 3;
to_col = num_of_layers + 2;
properties = {'flux', 'velocity_profile'};
%file_name = 'random_S1_11gy_11gx_15000WL_1_knick_9_output';
     [flux_data, header] = xlsread(strcat(dir, file_name,'.xlsx'), strcat('means_',properties{1}));
      [flux_std, header] = xlsread(strcat(dir, file_name,'.xlsx'), strcat('stds_',properties{1}));
%     [flux_stds, header] = xlsread(strcat(dir, file_name,'.xlsx'), strcat('stds_',properties{1}));
   %  [vel_data, header] = xlsread(strcat(dir, file_name,'.xlsx'), 'means_velocity_profile');
  %   [vel_stds, header] = xlsread(strcat(dir, file_name,'.xlsx'), 'stds_velocity_profile');
       % u(z=Hs) = 0
     Hs = h;%+h/(2*(num_of_layers));
   
    standa = strcat(num2str(6 * sqrt(g*solid_fraction)*(Hs)^(5/2)  / (15 *d)),'*(cosd(x)^0.5) *(tand(x) - a) /b');
    
    furbish = strcat(num2str(Hs),'*a*tand(x)/((1-tand(x)/b)*cosd(x))');
    caught_s = 0;
    caught_f = 0;
        %%%%  start flux  %%%%
    try
        [furb1_f, furb1_goodness] = fit(flux_data(1:transition_angle,1), flux_data(1:transition_angle, total_col), furbish);%,  'lower', [0, 0.1], 'upper', [10, 0.5]);
    catch e
        display('furb1_fit');
        display(file_name);
        e
        caught_f = 1;
        
    end
    while furb1_goodness.rsquare < 0.9
       % display(strcat('bad furbish fit-  R^2 = ', num2str(furb1_goodness.rsquare)));
       % furb1_f
        [furb1_f, furb1_goodness] = fit(flux_data(1:transition_angle,1), flux_data(1:transition_angle, total_col), furbish);%,  'lower', [0, 0.1], 'upper', [10, 0.5]);
    end
    
    try
       
             [st2_f, st2_goodness] = fit(flux_data(transition_angle:end,1), flux_data(transition_angle:end, total_col),standa);
       
     catch e
         display('st2_fit');
        display(file_name);
        e
        caught_s = 1;
    end
    if caught_f == 0
       
        K = furb1_f.a;
        Sc = furb1_f.b;
        fint = confint(furb1_f);
        stdK = std([fint(1,1), K, fint(2,1)]);
        stdSc = std([fint(1,2), Sc, fint(2,2)]);
         r_furb = furb1_goodness.rsquare;
    else
         K = -1;
        Sc = -1;
         r_furb = -1;
    end
    if caught_s == 0
   
        Mu = st2_f.a;
        Beta = st2_f.b;
        fint = confint(st2_f);
        stdMu = std([fint(1,1), Mu fint(2,1)]);
        stdBeta = std([fint(1,2), Beta,fint(2,2)]);
        r_standa = st2_goodness.rsquare;
    else
        Mu = -1;
        Beta = -1;
        r_standa = -1;
    end
    fluxes = flux_data(1:end, total_col);
    fluxes_stds = flux_std(1:end, total_col);
    Ek =  flux_data(1:end, 2);
    angles = flux_data(1:end, 1);
    trans_angle = num2str(flux_data(transition_angle,1));
    display(strcat(file_name, ', trans angle = ,',num2str(flux_data(transition_angle,1)), '  fu = ,', num2str(r_furb), ' st = ,', num2str(r_standa)));
end