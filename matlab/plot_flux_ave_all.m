M_pow = 2.5;
 WL_pow = 1;
colors = [0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8];
aves = zeros(2,16,13);
 aves(1,1:end,1:end) = xlsread('output/final/final_tables.xlsx', 'flux_S_average');
 aves(2,1:end,1:end) = xlsread('output/final/final_tables.xlsx', 'flux_std_S_average');
 aves(:,2:end,:) = aves(:,2:end,:) ;
 mwl035 = (aves(1,2:end, 1).^ M_pow) .* aves(1,2:end, 2) .^ WL_pow;
 figure('name', 'total flux');
 
 %subplot(1,2,1);
 hold on; 
 plots = zeros(2,11);
 for i = 3:13
     
     [fave, gave] = fit(mwl035.', aves(1,2:end,i).', 'a*x+b');
    % errorbar(mwl035.', aves(1,2:end,i).', aves(2,2:end,i).', '.', 'color', 'k');

    fave_plot = plot(fave, mwl035.', aves(1,2:end,i).');
    set(get(get(fave_plot(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(fave_plot(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    dplot = plot(mwl035.', aves(1,2:end,i).', 'o', 'MarkerSize',4, 'linewidth', 3);
    set(fave_plot(2), 'color', colors(i -2,1:3));
    set(dplot, 'color', colors(i -2,1:3), 'linewidth', 3);
    plots(:,i-2) = dplot;
    display(num2str(fave.b));
 end
 ylabel('Total Flux', 'fontsize', 14);
 xlabel(strcat('Magnitude^{', num2str(M_pow), '} * Duration'), 'fontsize', 14);
%legend off; 
 lh = legend(arrayfun(@num2str, aves(1,1,3:13), 'UniformOutput', false), 'location', 'southeast', 'FontSize',12);     

