plog=1;

colors = [0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8];

syses = {'S1', 'S2', 'S3'};
data = zeros(3,15,13);
aves = zeros(15,13);
mins = zeros(15,13);
maxs = zeros(15,13);
 [data(1,1:end,1:end), header] = xlsread('output/final/coefficients_std.xlsx', strcat('coeffs_',syses{1}));
 [data(2,1:end,1:end), header] = xlsread('output/final/coefficients_std.xlsx', strcat('coeffs_',syses{2}));
 [data(3,1:end,1:end), header] = xlsread('output/final/coefficients_std.xlsx', strcat('coeffs_',syses{3}));
  aves(:,:) = (data(1, :,:) + data(2, :,:) + data(3, :,:)) ./3;
  mins(:,:) = min(min(data(1, :,:),  data(2, :,:) ), data(3, :,:));
  maxs(:,:) = max(max(data(1, :,:),  data(2, :,:) ), data(3, :,:));
  unlog_aves = (data(1, :,:) + data(2, :,:) + data(3, :,:)) ./3;
  if plog == 1
      aves(:,:) = log(aves(:,:) );
      mins = aves-log(mins);
      maxs = log(maxs)-aves;
  else
      mins = aves-mins;
      maxs = maxs-aves;
  end
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
plot_by = K; 
plot_by_string='K';
if plot_by == 5
    plot_by_string='Sc';
end
if plot_by == 10
    plot_by_string='beta';
end
if plot_by == 8
    plot_by_string='mu';
end
WL_pow = 1;
M_pow = 1;
 % plot by WL
 figure('name', header{plot_by});
 
 %subplot(1,2,1);
 hold on;   
%plot data + std

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
%WL_pow = (log(aves(end-2,plot_by)) - log(aves(1,plot_by)))/(log(data(1,end-2,WL)) - log(data(1,1,WL)));

%     
%  e = errorbar(log(data(1,1:7,WL).^WL_pow), aves(1:7,plot_by), min(1:7,plot_by), max(1:7,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color','k');
%  e.LineStyle = 'none';
%      e = errorbar(log(data(1,8:10,WL).^WL_pow), aves(8:10,plot_by), min(8:10,plot_by),max(8:10,plot_by),'.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color', 'b');
%      e.LineStyle = 'none';
%      e = errorbar(log(data(1,11:13,WL).^WL_pow),aves(11:13,plot_by),min(11:13,plot_by),max(11:13,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1,'color', 'r');
%      e.LineStyle = 'none';
%      e = errorbar(log(data(1,14:15,WL).^WL_pow),aves(14:15,plot_by),min(14:15,plot_by),max(14:15,plot_by)'.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color', 'm');   
%      e.LineStyle = 'none';
     
     
if plog == 0 
     plot(data(1,1:7,WL).^WL_pow, aves(1:7,plot_by), 'color','k');
     plot(data(1,8:10,WL).^WL_pow, aves(8:10,plot_by), 'color', 'b');
     plot(data(1,11:13,WL).^WL_pow,aves(11:13,plot_by), 'color', 'r');
     plot(data(1,14:15,WL).^WL_pow,aves(14:15,plot_by), 'color', 'm');
     ylabel(header{plot_by});
     xlabel(strcat('Duration^{',num2str(WL_pow),'}'));
else
   

     
          [f1 g1] = fit(log(data(1,1:7,WL).^WL_pow).', aves(1:7,plot_by), 'a*x+b');
     [f2 g2] = fit(log(data(1,8:10,WL).^WL_pow).', aves(8:10,plot_by), 'a*x+b');
     [f3 g3] = fit(log(data(1,11:13,WL).^WL_pow).',aves(11:13,plot_by), 'a*x+b');
     [f4 g4] = fit(log(data(1,14:15,WL).^WL_pow).',aves(14:15,plot_by), 'a*x+b');
      f1_p = plot(f1,log(data(1,1:7,WL).^WL_pow).', aves(1:7,plot_by));
     f2_p = plot(f2,log(data(1,8:10,WL).^WL_pow).', aves(8:10,plot_by));
     f3_p = plot(f3,log(data(1,11:13,WL).^WL_pow).',aves(11:13,plot_by));
    % f4_p = plot(f4,log(data(1,14:15,WL).^WL_pow).',aves(14:15,plot_by));
   
 set(f1_p(2), 'color', 'k','LineWidth',3);
 set(f2_p(2), 'color', 'b','LineWidth',3);
 set(f3_p(2), 'color', 'r','LineWidth',3);
 %set(f4_p(2), 'color', 'm','LineWidth',3);
 
      plot(log(data(1,1:7,WL).^WL_pow), aves(1:7,plot_by), 'o','MarkerSize',3,'MarkerEdgeColor','k', 'color','k');
     plot(log(data(1,8:10,WL).^WL_pow), aves(8:10,plot_by),'o','MarkerSize',3,'MarkerEdgeColor','b', 'color', 'b');
     plot(log(data(1,11:13,WL).^WL_pow),aves(11:13,plot_by), 'o','MarkerSize',3,'MarkerEdgeColor','r','color', 'r');
    % plot(log(data(1,14:15,WL).^WL_pow),aves(14:15,plot_by),'.', 'color', 'm');
 
   e = errorbar(log(data(1,1:7,WL).^WL_pow), aves(1:7,plot_by), mins(1:7,plot_by), maxs(1:7,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color','k');
%  e.LineStyle = 'none';
      e = errorbar(log(data(1,8:10,WL).^WL_pow), aves(8:10,plot_by), mins(8:10,plot_by),maxs(8:10,plot_by),'.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color', 'b');
%      e.LineStyle = 'none';
      e = errorbar(log(data(1,11:13,WL).^WL_pow),aves(11:13,plot_by),mins(11:13,plot_by),maxs(11:13,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1,'color', 'r');
%      e.LineStyle = 'none';
%      e = errorbar(log(data(1,14:15,WL).^WL_pow),aves(14:15,plot_by),min(14:15,plot_by),max(14:15,plot_by)'.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color', 'm');   
%      e.LineStyle = 'none';

% text(8.5, -5.5,...
%     {strcat(num2str(f1.a),'X + ',num2str(f1.b), ' , R^2 = ', num2str(g1.rsquare)),...
%     strcat(num2str(f2.a),'X + ',num2str(f2.b), ' , R^2 = ', num2str(g2.rsquare)),...
%     strcat(num2str(f3.a),'X + ',num2str(f3.b), ' , R^2 = ', num2str(g3.rsquare)),...
%     strcat(num2str(f4.a),'X + ',num2str(f4.b), ' , R^2 = ', num2str(g4.rsquare))},...
%     'HorizontalAlignment','left',...
%           'VerticalAlignment','top',...
%           'BackgroundColor',[1 1 1],...
%           'FontSize',14);
     ylabel(strcat('log   \', header{plot_by}), 'FontSize',14);
     xlabel(strcat('log Duration'), 'FontSize',14);   
  for i = (14:0.5:17)
   %   for i = (-7.3:-0.05:-7.38)
 plot(log(data(1,1:13,WL)),log(data(1,1:13,WL))-i, ':');
    end
 legend([f1_p(2) f2_p(2) f3_p(2)],strcat('mag = ', num2str(data(1,7,Mag)),'    R^2= ',num2str(g1.rsquare)),...
     strcat('mag = ', num2str(data(1,10,Mag)),'    R^2= ',num2str(g2.rsquare)),...
     strcat('mag = ', num2str(data(1,13,Mag)),'    R^2= ',num2str(g3.rsquare)),...
     'font', 14, 'location', 'northwest');
  
     %legend({'Magnitude=11    ', num2str(g1.rsquare); 'Magnitude=19    ' ,num2str(g2.rsquare); 'Magnitude=25    ', num2str(g3.rsquare)}, 'location', 'northwest');  
end
 
  
     display(strcat(num2str(f1.a), '  mag = ', num2str(data(1,7,Mag)),' , ', num2str((log(aves(7,plot_by)) - log(aves(1,plot_by)))/(log(data(1,7,WL)) - log(data(1,1,WL))))));
     display(strcat(num2str(f2.a), '  mag = ', num2str(data(1,10,Mag)),' , ', num2str((log(aves(10,plot_by)) - log(aves(8,plot_by)))/(log(data(1,10,WL)) - log(data(1,8,WL))))));
     display(strcat(num2str(f3.a), '  mag = ', num2str(data(1,13,Mag)),' , ', num2str((log(aves(13,plot_by)) - log(aves(11,plot_by)))/(log(data(1,13,WL)) - log(data(1,11,WL))))));
     display(strcat('mag = ', num2str(data(1,15,Mag)),' , ', num2str((log(aves(15,plot_by)) - log(aves(14,plot_by)))/(log(data(1,15,WL)) - log(data(1,14,WL))))));
     
     
     saveas(gcf, strcat('output/final2/ave_',plot_by_string,'_std_by_mag.fig'));
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
  aves = sortrows(aves, [2 1]);
 
  %subplot(1,2,2);
  figure('name', header{plot_by});
    hold on;
     
%      errorbar(data(1,3:5,Mag),data(1,3:5,K),data(1,3:5,stdK), 'color', colors(1, 1:3));
%      errorbar(data(1,8:10,Mag),data(1,8:10,K),data(1,8:10,stdK), 'color', colors(2, 1:3));
%      errorbar(data(1,11:13,Mag),data(1,11:13,K),data(1,11:13,stdK), 'color', colors(3, 1:3));
%      errorbar(data(2,3:5,Mag),data(2,3:5,K),data(2,3:5,stdK), 'color', colors(4, 1:3));
%      errorbar(data(2,8:10,Mag),data(2,8:10,K),data(2,8:10,stdK), 'color', colors(5, 1:3));
%      errorbar(data(2,11:13,Mag),data(2,11:13,K),data(2,11:13,stdK), 'color', colors(6, 1:3));
%      errorbar(data(2,3:5,Mag),data(2,3:5,K),data(2,3:5,stdK), 'color', colors(7, 1:3));
%      errorbar(data(2,8:10,Mag),data(2,8:10,K),data(2,8:10,stdK), 'color', colors(8, 1:3));
%      errorbar(data(2,11:13,Mag),data(2,11:13,K),data(2,11:13,stdK), 'color', colors(9, 1:3));
%M_pow = (log(aves(end-2,plot_by)) - log(aves(1,plot_by)))/(log(s1_mag(1,end-2,Mag)) - log(s1_mag(1,1,Mag)));

if plog == 0
%      plot(s1_mag(3:6,Mag).^M_pow,aves(3:6,plot_by), 'color','k');
%      plot(s1_mag(7:10,Mag).^M_pow,aves(7:10,plot_by), 'color','b');
%      plot(s1_mag(11:13,Mag).^M_pow,aves(11:13,plot_by), 'color', 'r');
          [f1 g1] = fit(s1_mag(3:6,Mag), aves(3:6,plot_by), 'a*x+b');
     [f2 g2] = fit(s1_mag(7:10,Mag), aves(7:10,plot_by), 'a*x+b');
     [f3 g3] = fit(s1_mag(11:13,Mag),aves(11:13,plot_by), 'a*x+b');
      f1_p = plot(f1,s1_mag(3:6,Mag), aves(3:6,plot_by));
     f2_p = plot(f2,s1_mag(7:10,Mag), aves(7:10,plot_by));
     f3_p = plot(f3,s1_mag(11:13,Mag),aves(11:13,plot_by));
   
 set(f1_p(2), 'color', 'k','LineWidth',3);
 set(f2_p(2), 'color', 'b','LineWidth',3);
 set(f3_p(2), 'color', 'r','LineWidth',3);
      e = errorbar(s1_mag(3:6,Mag), aves(3:6,plot_by), mins(3:6,plot_by), maxs(3:6,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color','k');
%  e.LineStyle = 'none';
      e = errorbar(s1_mag(7:10,Mag), aves(7:10,plot_by), mins(7:10,plot_by),maxs(7:10,plot_by),'.','MarkerSize',3,'MarkerEdgeColor','b','LineWidth',1, 'color', 'b');
%      e.LineStyle = 'none';
      e = errorbar(s1_mag(11:13,Mag),aves(11:13,plot_by),mins(11:13,plot_by),maxs(11:13,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','r','LineWidth',1,'color', 'r');
     ylabel('\mu_d (non-dimensional)','FontSize',14);
     xlabel(strcat('Magnitude (non-dimensional)'),'FontSize',14);
      legend([f1_p(2) f2_p(2) f3_p(2)],...
         strcat('Duration = ', num2str(s1_mag(3,WL)),' slope = ',num2str(f1.a),'  R^2= ',num2str(g1.rsquare)),...
     strcat('dur = ', num2str(s1_mag(7,WL)),' slope = ',num2str(f2.a),'  R^2= ',num2str(g2.rsquare)),...
     strcat('dur = ', num2str(s1_mag(11,WL)),' slope = ',num2str(f3.a),'  R^2= ',num2str(g3.rsquare)),...
     'location', 'northeast');
else
     
     
          [f1 g1] = fit(log(s1_mag(3:6,Mag).^M_pow), aves(3:6,plot_by), 'a*x+b');
     [f2 g2] = fit(log(s1_mag(7:10,Mag).^M_pow), aves(7:10,plot_by), 'a*x+b');
     [f3 g3] = fit(log(s1_mag(11:13,Mag).^M_pow),aves(11:13,plot_by), 'a*x+b');
      f1_p = plot(f1,log(s1_mag(3:6,Mag).^M_pow), aves(3:6,plot_by));
     f2_p = plot(f2,log(s1_mag(7:10,Mag).^M_pow), aves(7:10,plot_by));
     f3_p = plot(f3,log(s1_mag(11:13,Mag).^M_pow),aves(11:13,plot_by));
   
 set(f1_p(2), 'color', 'k','LineWidth',3);
 set(f2_p(2), 'color', 'b','LineWidth',3);
 set(f3_p(2), 'color', 'r','LineWidth',3);
%     plot(log(s1_mag(3:6,Mag).^M_pow),aves(3:6,plot_by), 'o','color','k');
%      plot(log(s1_mag(7:10,Mag).^M_pow),aves(7:10,plot_by),'o', 'color','b');
%      plot(log(s1_mag(11:13,Mag).^M_pow),aves(11:13,plot_by),'o', 'color', 'r');
    %for i = (12:0.5:14)
    for i = (12.5:0.3:13.5)
 plot(log(s1_mag(3:13,Mag)),2.5*log(s1_mag(3:13,Mag) )-i, ':');
 %plot(log(s1_mag(3:13,Mag).^M_pow),log(s1_mag(3:13,Mag) .^M_pow)-i, ':');
    end
   legend([f1_p(2) f2_p(2) f3_p(2)],...
         strcat('Duration = ', num2str(s1_mag(3,WL)),' slope = ',num2str(f1.a),'  R^2= ',num2str(g1.rsquare)),...
     strcat('Duration = ', num2str(s1_mag(7,WL)),' slope = ',num2str(f2.a),'  R^2= ',num2str(g2.rsquare)),...
     strcat('Duration = ', num2str(s1_mag(11,WL)),' slope = ',num2str(f3.a),'  R^2= ',num2str(g3.rsquare)),...
     'location', 'northeast');
 
 
%      text(2.2, -5.5,...
%     {strcat(num2str(f1.a),'X + ',num2str(f1.b), ' , R^2 = ', num2str(g1.rsquare)),...
%     strcat(num2str(f2.a),'X + ',num2str(f2.b), ' , R^2 = ', num2str(g2.rsquare)),...
%     strcat(num2str(f3.a),'X + ',num2str(f3.b), ' , R^2 = ', num2str(g3.rsquare))},.....
%     'HorizontalAlignment','left',...
%           'VerticalAlignment','top',...
%           'BackgroundColor',[1 1 1],...
%           'FontSize',14);
     e = errorbar(log(s1_mag(3:6,Mag).^M_pow), aves(3:6,plot_by), mins(3:6,plot_by), maxs(3:6,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','k','LineWidth',1, 'color','k');
%  e.LineStyle = 'none';
      e = errorbar(log(s1_mag(7:10,Mag).^M_pow), aves(7:10,plot_by), mins(7:10,plot_by),maxs(7:10,plot_by),'.','MarkerSize',3,'MarkerEdgeColor','b','LineWidth',1, 'color', 'b');
%      e.LineStyle = 'none';
      e = errorbar(log(s1_mag(11:13,Mag).^M_pow),aves(11:13,plot_by),mins(11:13,plot_by),maxs(11:13,plot_by), '.','MarkerSize',3,'MarkerEdgeColor','r','LineWidth',1,'color', 'r');
    ylabel(strcat('log  ', header{plot_by}), 'FontSize',14);
     xlabel(strcat('log Magnitude'),'FontSize', 14);%^{', num2str(M_pow),'}'), 'FontSize',14);
end
%inferred best fit slope'
%     legend({'Duration=10000    ' num2str(g1.rsquare); 'Duration=12000    ' num2str(g2.rsquare); 'Duration=15000    ' num2str(g3.rsquare)}, 'location', 'northwest');

     display(strcat(num2str(f1.a), '  wl = ', num2str(s1_mag(6,WL)),' , ', num2str((log(aves(6,plot_by)) - log(aves(3,plot_by)))/(log(s1_mag(6,Mag)) - log(s1_mag(3,Mag))))));
     display(strcat(num2str(f2.a), '  wl = ', num2str(s1_mag(10,WL)),' , ', num2str((log(aves(10,plot_by)) - log(aves(7,plot_by)))/(log(s1_mag(10,Mag)) - log(s1_mag(7,Mag))))));
     display(strcat(num2str(f3.a), '  wl = ', num2str(s1_mag(13,WL)),' , ', num2str((log(aves(13,plot_by)) - log(aves(11,plot_by)))/(log(s1_mag(13,Mag)) - log(s1_mag(11,Mag))))));
     
     saveas(gcf, strcat('output/final2/ave_',plot_by_string,'_std_by_wl.fig'));
 %plot by f(mag,WL)   
 