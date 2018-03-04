warning('off','all');
colors = [0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8];

syses = {'S1', 'S2', 'S3'};
data = zeros(3, 16,13);
stds = zeros(3, 16,13);
cols = [13 12 13];
 [data(1,1:end,1:cols(1))] = xlsread('output/final/coefficients_std.xlsx', strcat('flux_',syses{1}));
 [data(2,1:end,1:cols(2))] = xlsread('output/final/coefficients_std.xlsx', strcat('flux_',syses{2}));
 [data(3,1:end,1:cols(3))] = xlsread('output/final/coefficients_std.xlsx', strcat('flux_',syses{3}));
 [stds(1,1:end,1:cols(1))] = xlsread('output/final/coefficients_std.xlsx', strcat('flux_std_',syses{1}));
 [stds(2,1:end,1:cols(2))] = xlsread('output/final/coefficients_std.xlsx', strcat('flux_std_',syses{2}));
 [stds(3,1:end,1:cols(3))] = xlsread('output/final/coefficients_std.xlsx', strcat('flux_std_',syses{3}));
 M_pow = 2.5;
 WL_pow = 1;
  Mag=1;
 WL=2;
mwl035 = (data(1,2:end, Mag).^ M_pow);% .* data(1,2:end, WL) .^ WL_pow;
 % plot by WL
 for i = 1:numel(syses)
     figure('name', strcat('flux vs MWL ', syses{i}));
     hold on;  
     for j = 3:cols(i)-1
       f = fit(mwl035.', data(i,2:end,j).', 'a*x+b');
       f_plot = plot(f, mwl035, data(i,2:end,j));
        set(f_plot(2), 'color', colors(j, 1:3));
         errorbar(mwl035.', data(i,2:end,j).', stds(i,2:end,j), '.', 'color', colors(j, 1:3));
        % plot(mwl035.', data(3,1:end,end).',  '.', 'color', 'b');
     end

     ylabel('flux');
      xlabel(strcat('M^{',num2str(M_pow),'}*WL^{', num2str(WL_pow),'}'));
      legend off;
 end
 