
colors = [0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8];

syses = {'S1', 'S2', 'S3'};
data = zeros(3,15,13);
 [data(1,1:end,1:end), header] = xlsread('output/final/coefficients_std.xlsx', strcat('coeffs_',syses{1}));
 [data(2,1:end,1:end), header] = xlsread('output/final/coefficients_std.xlsx', strcat('coeffs_',syses{2}));
 [data(3,1:end,1:end), header] = xlsread('output/final/coefficients_std.xlsx', strcat('coeffs_',syses{3}));
 
 Mag=1;
 WL=2;
 K=3;
 stdK =4;
 Sc =5;
 stdSc=6;
 fR2= 7;
 mu_d = 8;
 stdMu = 9;
 beta = 10;
 stdBeta= 11;
 sR2 = 12;
 trangle = 13;
plot_by = mu_d; 
WL_pow = 1;
M_pow = 1;
 % plot by WL
 figure('name', header{plot_by});
 
 subplot(1,2,1);
 hold on;   
 ylabel(header{plot_by});
 xlabel(strcat('WL^{',num2str(WL_pow),'}'));
%      errorbar(data(1,1:7,WL),data(1,1:7,K),data(1,1:7,stdK), 'color',colors(1, 1:3));
%      errorbar(data(1,8:10,WL),data(1,8:10,K),data(1,8:10,stdK), 'color', colors(2, 1:3));
%      errorbar(data(1,11:13,WL),data(1,11:13,K),data(1,11:13,stdK), 'color', colors(3, 1:3));
%      errorbar(data(2,1:7,WL),data(2,1:7,K),data(2,1:7,stdK), 'color', colors(4, 1:3));
%      errorbar(data(2,8:10,WL),data(2,8:10,K),data(2,8:10,stdK), 'color', colors(5, 1:3));
%      errorbar(data(2,11:13,WL),data(2,11:13,K),data(2,11:13,stdK), 'color', colors(6, 1:3));
%      errorbar(data(3,1:7,WL),data(3,1:7,K),data(3,1:7,stdK), 'color', colors(7, 1:3));
%      errorbar(data(3,8:10,WL),data(3,8:10,K),data(3,8:10,stdK), 'color', colors(8, 1:3));
%      errorbar(data(3,11:13,WL),data(3,11:13,K),data(3,11:13,stdK), 'color', colors(9, 1:3));
% % 
     plot(data(1,1:7,WL).^WL_pow,data(1,1:7,plot_by), 'color',colors(1, 1:3));
     plot(data(1,8:10,WL).^WL_pow,data(1,8:10,plot_by), 'color', colors(2, 1:3));
     plot(data(1,11:13,WL).^WL_pow,data(1,11:13,plot_by), 'color', colors(3, 1:3));
     plot(data(1,14:15,WL).^WL_pow,data(1,14:15,plot_by), 'color', colors(4, 1:3));
     plot(data(2,1:7,WL).^WL_pow,data(2,1:7,plot_by), 'color', colors(5, 1:3));
     plot(data(2,8:10,WL).^WL_pow,data(2,8:10,plot_by), 'color', colors(6, 1:3));
     plot(data(2,11:13,WL).^WL_pow,data(2,11:13,plot_by), 'color', colors(7, 1:3));
     plot(data(2,14:15,WL).^WL_pow,data(2,14:15,plot_by), 'color', colors(8, 1:3));
     plot(data(3,1:7,WL).^WL_pow,data(3,1:7,plot_by), 'color', colors(9, 1:3));
     plot(data(3,8:10,WL).^WL_pow,data(3,8:10,plot_by), 'color', colors(10, 1:3));
     plot(data(3,11:13,WL).^WL_pow,data(3,11:13,plot_by), 'color', colors(11, 1:3));
     plot(data(3,14:15,WL).^WL_pow,data(3,14:15,plot_by), 'color', colors(12, 1:3));

 legend([strcat(syses{1}, '-', num2str(data(1,3,Mag))) ; strcat(syses{1}, '-', num2str(data(1,8,Mag))) ; strcat(syses{1}, '-', num2str(data(1,11,Mag))) ; strcat(syses{1}, '-', num2str(data(1,14,Mag))) ; ...
            strcat(syses{2}, '-', num2str(data(1,3,Mag))) ; strcat(syses{2}, '-', num2str(data(1,8,Mag))) ; strcat(syses{2}, '-', num2str(data(1,11,Mag))) ; strcat(syses{2}, '-', num2str(data(1,14,Mag))) ;...
            strcat(syses{3}, '-', num2str(data(1,3,Mag))) ; strcat(syses{3}, '-', num2str(data(1,8,Mag))) ; strcat(syses{3}, '-', num2str(data(1,11,Mag))) ; strcat(syses{3}, '-', num2str(data(1,14,Mag))) ;] , ...
            'location', 'northeast');
  
 % plot by Mag
 s1_mag = zeros(15,13);
 s2_mag = zeros(15,13);
 s3_mag = zeros(15,13);
 s1_mag(:,:) = data(1,1:end,1:end);
 s2_mag(:,:) =data(2,1:end,1:end);
 s3_mag(:,:) =data(3,1:end,1:end);
 s1_mag = sortrows(s1_mag, [2 1]);
 s2_mag = sortrows(s2_mag, [2 1]);
 s3_mag = sortrows(s3_mag, [2 1]);
 
  subplot(1,2,2);
    hold on;
     ylabel(header{plot_by});
     xlabel(strcat('Mag^{', num2str(M_pow),'}'));
%      errorbar(data(1,3:5,Mag),data(1,3:5,K),data(1,3:5,stdK), 'color', colors(1, 1:3));
%      errorbar(data(1,8:10,Mag),data(1,8:10,K),data(1,8:10,stdK), 'color', colors(2, 1:3));
%      errorbar(data(1,11:13,Mag),data(1,11:13,K),data(1,11:13,stdK), 'color', colors(3, 1:3));
%      errorbar(data(2,3:5,Mag),data(2,3:5,K),data(2,3:5,stdK), 'color', colors(4, 1:3));
%      errorbar(data(2,8:10,Mag),data(2,8:10,K),data(2,8:10,stdK), 'color', colors(5, 1:3));
%      errorbar(data(2,11:13,Mag),data(2,11:13,K),data(2,11:13,stdK), 'color', colors(6, 1:3));
%      errorbar(data(2,3:5,Mag),data(2,3:5,K),data(2,3:5,stdK), 'color', colors(7, 1:3));
%      errorbar(data(2,8:10,Mag),data(2,8:10,K),data(2,8:10,stdK), 'color', colors(8, 1:3));
%      errorbar(data(2,11:13,Mag),data(2,11:13,K),data(2,11:13,stdK), 'color', colors(9, 1:3));

     plot(s1_mag(3:6,Mag).^M_pow,s1_mag(3:6,plot_by), 'color', colors(1, 1:3));
     plot(s1_mag(7:10,Mag).^M_pow,s1_mag(7:10,plot_by), 'color', colors(2, 1:3));
     plot(s1_mag(11:13,Mag).^M_pow,s1_mag(11:13,plot_by), 'color', colors(3, 1:3));
     plot(s2_mag(3:6,Mag).^M_pow,s2_mag(3:6,plot_by), 'color', colors(4, 1:3));
     plot(s2_mag(7:10,Mag).^M_pow,s2_mag(7:10,plot_by), 'color', colors(5, 1:3));
     plot(s2_mag(11:13,Mag).^M_pow,s2_mag(11:13,plot_by), 'color', colors(6, 1:3));
     plot(s3_mag(3:6,Mag).^M_pow,s3_mag(3:6,plot_by), 'color', colors(7, 1:3));
     plot(s3_mag(7:10,Mag).^M_pow,s3_mag(7:10,plot_by), 'color', colors(8, 1:3));
     plot(s3_mag(11:13,Mag).^M_pow,s3_mag(11:13,plot_by), 'color', colors(9, 1:3));
 
 legend([strcat(syses{1}, '-', num2str(s1_mag(3,WL))) ; strcat(syses{1}, '-', num2str(s1_mag(7,WL))) ; strcat(syses{1}, '-', num2str(s1_mag(11,WL))) ; ...
            strcat(syses{2}, '-', num2str(s1_mag(3,WL))) ; strcat(syses{2}, '-', num2str(s1_mag(7,WL))) ; strcat(syses{2}, '-', num2str(s1_mag(11,WL))) ; ...
            strcat(syses{3}, '-', num2str(s1_mag(3,WL))) ; strcat(syses{3}, '-', num2str(s1_mag(7,WL))) ; strcat(syses{3}, '-', num2str(s1_mag(11,WL))) ;], ...
            'location', 'southeast');
  
 %plot by f(mag,WL)   
 