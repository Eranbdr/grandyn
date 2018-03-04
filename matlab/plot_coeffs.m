
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
 tangle = 13;
 plots = [K, Sc, mu_d, beta];
 
 
 % plot by WL
 figure('name', 'by WL');
 title('WL');
 
 for i = plots(1:end)
     subplot(2,2,find(plots==i,1));
     hold on;   
     ylabel(header(i));
     xlabel('WL');
%      errorbar(data(1,1:7,WL),data(1,1:7,i),data(1,1:7,i+1), '.', 'color',colors(1, 1:3));
%      errorbar(data(1,8:10,WL),data(1,8:10,i),data(1,8:10,i+1), '.', 'color', colors(2, 1:3));
%      errorbar(data(1,11:13,WL),data(1,11:13,i),data(1,11:13,i+1), '.', 'color', colors(3, 1:3));
%      errorbar(data(2,1:7,WL),data(2,1:7,i),data(2,1:7,i+1), '.', 'color', colors(4, 1:3));
%      errorbar(data(2,8:10,WL),data(2,8:10,i),data(2,8:10,i+1), '.', 'color', colors(5, 1:3));
%      errorbar(data(2,11:13,WL),data(2,11:13,i),data(2,11:13,i+1), '.', 'color', colors(6, 1:3));

     plot(data(1,1:7,WL),data(1,1:7,i), '.', 'color',colors(1, 1:3));
     plot(data(1,8:10,WL),data(1,8:10,i), '.', 'color', colors(2, 1:3));
     plot(data(1,11:13,WL),data(1,11:13,i), '.', 'color', colors(3, 1:3));
     plot(data(1,14:15,WL),data(1,14:15,i), '.', 'color', colors(4, 1:3));
     plot(data(2,1:7,WL),data(2,1:7,i), '.', 'color', colors(5, 1:3));
     plot(data(2,8:10,WL),data(2,8:10,i), '.', 'color', colors(6, 1:3));
     plot(data(2,11:13,WL),data(2,11:13,i), '.', 'color', colors(7, 1:3));
     plot(data(2,14:15,WL),data(2,14:15,i), '.', 'color', colors(8, 1:3));
     plot(data(3,1:7,WL),data(3,1:7,i), '.', 'color', colors(9, 1:3));
     plot(data(3,8:10,WL),data(3,8:10,i), '.', 'color', colors(10, 1:3));
     plot(data(3,11:13,WL),data(3,11:13,i), '.', 'color', colors(11, 1:3));
     plot(data(3,14:15,WL),data(3,14:15,i), '.', 'color', colors(12, 1:3));
 end
 legend([strcat(syses{1}, '-', num2str(data(1,3,Mag))) ; strcat(syses{1}, '-', num2str(data(1,8,Mag))) ; strcat(syses{1}, '-', num2str(data(1,11,Mag))) ; strcat(syses{1}, '-', num2str(data(1,14,Mag))) ; ...
            strcat(syses{2}, '-', num2str(data(1,3,Mag))) ; strcat(syses{2}, '-', num2str(data(1,8,Mag))) ; strcat(syses{2}, '-', num2str(data(1,11,Mag))) ; strcat(syses{2}, '-', num2str(data(1,14,Mag))) ;...
            strcat(syses{3}, '-', num2str(data(1,3,Mag))) ; strcat(syses{3}, '-', num2str(data(1,8,Mag))) ; strcat(syses{3}, '-', num2str(data(1,11,Mag))) ; strcat(syses{3}, '-', num2str(data(1,14,Mag))) ;] , ...
            'location', 'southeast');
  
 % plot by Mag
 figure('name', 'by Mag');
 title('Mag');
 s1_mag = zeros(15,13);
 s2_mag = zeros(15,13);
 s3_mag = zeros(15,13);
 s1_mag(:,:) = data(1,1:end,1:end);
 s2_mag(:,:) =data(2,1:end,1:end);
 s3_mag(:,:) =data(3,1:end,1:end);
 s1_mag = sortrows(s1_mag, [2 1]);
 s2_mag = sortrows(s2_mag, [2 1]);
 s3_mag = sortrows(s3_mag, [2 1]);
 
 for i = plots(1:end)
     subplot(2,2,find(plots==i,1));
     hold on;   
     ylabel(header(i));
     xlabel('Mag');
%      errorbar(data(1,3:5,Mag),data(1,3:5,i),data(1,3:5,i+1), '.', 'color', colors(1, 1:3));
%      errorbar(data(1,8:10,Mag),data(1,8:10,i),data(1,8:10,i+1), '.', 'color', colors(2, 1:3));
%      errorbar(data(1,11:13,Mag),data(1,11:13,i),data(1,11:13,i+1), '.', 'color', colors(3, 1:3));
%      errorbar(data(2,3:5,Mag),data(2,3:5,i),data(2,3:5,i+1), '.', 'color', colors(4, 1:3));
%      errorbar(data(2,8:10,Mag),data(2,8:10,i),data(2,8:10,i+1), '.', 'color', colors(5, 1:3));
%      errorbar(data(2,11:13,Mag),data(2,11:13,i),data(2,11:13,i+1), '.', 'color', colors(6, 1:3));

     plot(s1_mag(3:6,Mag),s1_mag(3:6,i),'.', 'color', colors(1, 1:3));
     plot(s1_mag(7:10,Mag),s1_mag(7:10,i), '.', 'color', colors(2, 1:3));
     plot(s1_mag(11:13,Mag),s1_mag(11:13,i), '.', 'color', colors(3, 1:3));
     plot(s2_mag(3:6,Mag),s2_mag(3:6,i), '.', 'color', colors(4, 1:3));
     plot(s2_mag(7:10,Mag),s2_mag(7:10,i), '.', 'color', colors(5, 1:3));
     plot(s2_mag(11:13,Mag),s2_mag(11:13,i), '.', 'color', colors(6, 1:3));
     plot(s3_mag(3:6,Mag),s3_mag(3:6,i), '.', 'color', colors(4, 1:3));
     plot(s3_mag(7:10,Mag),s3_mag(7:10,i), '.', 'color', colors(5, 1:3));
     plot(s3_mag(11:13,Mag),s3_mag(11:13,i), '.', 'color', colors(6, 1:3));
 end
 legend([strcat(syses{1}, '-', num2str(s1_mag(3,WL))) ; strcat(syses{1}, '-', num2str(s1_mag(7,WL))) ; strcat(syses{1}, '-', num2str(s1_mag(11,WL))) ; ...
            strcat(syses{2}, '-', num2str(s1_mag(3,WL))) ; strcat(syses{2}, '-', num2str(s1_mag(7,WL))) ; strcat(syses{2}, '-', num2str(s1_mag(11,WL))) ; ...
            strcat(syses{3}, '-', num2str(s1_mag(3,WL))) ; strcat(syses{3}, '-', num2str(s1_mag(7,WL))) ; strcat(syses{3}, '-', num2str(s1_mag(11,WL))) ;], ...
            'location', 'southeast');
  
 %plot by f(mag,WL)   
 