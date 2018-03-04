% read several CSVs and present means and vars 
% input is slope, file, SS_step
warning('off','all');
solid_fraction = 0.8;%965/(34*34);
g = 9.81;
d = 0.001;
rho_p = 2500;
rho = rho_p * solid_fraction;
depth = 34*0.001;
sub_name = 'random_S2_19g_12000WL_4';
columns = 13;% S 13, M 23
knick = 9; % from lowest slope
set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
colors = get(0, 'DefaultAxesColorOrder');

properties = {'flux', 'velocity_profile'};

[num, files] = xlsread(strcat('C:\Users\EranBD\Documents\thesis\matlab\', sub_name, '_input.csv'));
output_file = strcat('C:\Users\EranBD\Documents\thesis\matlab_2\output\', sub_name, '_knick_', num2str(knick),  '_output_1_2d_.xlsx');
delete output_file;
files_num = numel(files);
figure('name', strcat(sub_name, '_knick_', num2str(knick)));


    
    [header, means_t, vars_t, stds_t] = single_csv_data(files{1}, num(1,3), num(1,4), properties{1}, 1, columns);
    cols = numel(header);
    means= zeros(files_num, cols);
    vars = zeros(files_num, cols);
    stds = zeros(files_num, cols);
    means(1,1:end) = means_t;
    vars(1,1:end) = vars_t;
    stds(1,1:end) = stds_t;
    
    for d = 2:files_num
        [header, means(d,1:end), vars(d,1:end), stds(d,1:end)] = single_csv_data(files{d}, num(d,3), num(d,4),properties{1}, 1, columns);
    end
     means(1:end,1)= num(1:end,1);
     vars(1:end,1)= num(1:end,1);
     stds(1:end,1)= num(1:end,1);

     col = cols;
        

        %axis([0,0,tand(16),10]) ;
    
        subplot(2,1,1);
hold on;        
        [roer_f, roer_goodness] = fit(means(1:end,1), means(1:end, col),'(a*tand(x))/(1-(tand(x)/b)^2)','lower', [0,0.26]', 'upper', [500,0.4]);
       % [exp_f, exp_goodness] = fit(means(1:end,1), means(1:end, col),'exp1');
        roer_plot = plot(roer_f, means(1:end,1), means(1:end, col));
       % exp_plot = plot(exp_f, means(1:end,1), means(1:end, col));
        
        set(roer_plot(2), 'color', 'r');
        if knick > 2
            [roer1_f, roer1_goodness] = fit(means(1:knick,1), means(1:knick, col),'(a*tand(x))/(1-(tand(x)/b)^2)', 'lower', [0, 0.26], 'upper', [500, 0.3]);
             roer1_plot = plot(roer1_f, means(1:knick,1), means(1:knick, col));

          %  [exp1_f, exp1_goodness] = fit(means(1:r_knick,1), means(1:r_knick, col),'exp1');
            %exp1_plot = plot(exp1_f, means(1:r_knick,1), means(1:r_knick, col));
            
            set(roer1_plot(2), 'color', 'k');
        else
            [roer1_f, roer1_goodness] = fit((1:10).', (1:10).','(a*tand(x))/(1-(tand(x)/b)^2)', 'lower', [0, 0.26], 'upper', [500, 0.3]);
            %,'startpoint', [1e-7,means(knick,1)*0.9999], 'lower', [0,means(knick,1)*0.9], 'upper', [1, means(knick,1)*1.1]);
          
            [exp1_f, exp1_goodness] = fit((1:10).', (1:10).','exp1');
        end
        
        if knick < numel(files)
           
            standa = strcat(num2str(2 * rho * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_d) * ', num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *b)');
            [st2_f, st2_goodness] = fit(means(knick:end,1), means(knick:end, col),standa);
            st2_plot = plot(st2_f, means(knick:end,1), means(knick:end, col));
            
             %fernandez = strcat(num2str(2 * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_s) * I * ', num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *(m_2  - tand(x)))');
            %[fern2_f, fern2_goodness] = fit(means(knick:end,1), means(knick:end, col),fernandez, 'lower', [0 0.2 0], 'upper', [100 0.5 0.4]);
            %fern2_plot = plot(fern2_f, means(knick:end,1), means(knick:end, col));

            
            %set(fern2_plot(2), 'color', 'b');
            set(st2_plot(2), 'color', 'b');
        else
            %fernandez = strcat(num2str(2 * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_s) * I * ', num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *(m_2 + m_s - tand(x)))');
            %[fern2_f, fern2_goodness] = fit((1:10).',(1:10).',fernandez);
    
            standa = strcat(num2str(2 * rho * sqrt(solid_fraction*g)),'*(cosd(x)^0.5) *(tand(x) - m_d) * ', num2str(3/5*(34*0.001)^(5/2)),'  / (3 *0.001 *b)');
            [st2_f, st2_goodness] = fit((1:10).',(1:10).',standa);
        end
        
        
        
     e = errorbar(means(1:end,1),means(1:end, col),stds(1:end, col), '.'); 
    e.LineStyle = 'none';
        title('Total Flux');
        xlabel('slope angle');  
        ylabel('average total flux');
        legend('location', 'northwest');
        %strcat('Fernandez (b): m_s=',num2str(fern2_f.m_s), ' m_2=',num2str(fern2_f.m_2), ' I_0=', num2str(fern2_f.I),' rsquare=', num2str(fern2_goodness.rsquare))},...
        % strcat('Exp (k): A=',num2str(exp1_f.a), ' EXP=', num2str(exp1_f.b),' rsquare=', num2str(exp1_goodness.rsquare)),...
        %strcat('Standa (b): m_d=',num2str(st2_f.b), ' beta=', num2str(st2_f.b),' rsquare=', num2str(st2_goodness.rsquare)),...
        text((means(1,1)+means(floor(files_num/2),1))/4,means(end, col),...
          {strcat('Low Roering (k): K=',num2str(roer1_f.a), ' Sc=', num2str(atand(roer1_f.b)),' rsquare=', num2str(roer1_goodness.rsquare)),...                 
                    strcat('Standa (b): m_d=',num2str(st2_f.m_d), ' beta=', num2str(st2_f.b),' rsquare=', num2str(st2_goodness.rsquare)),...
                    strcat('Roering (r): K=',num2str(roer_f.a), ' Sc=', num2str(atand(roer_f.b)),' rsquare=', num2str(roer_goodness.rsquare))},...
          'HorizontalAlignment','left',...
          'VerticalAlignment','top',...
          'BackgroundColor',[1 1 1],...
          'FontSize',10);
     
    legend('off');
      xlswrite(output_file, header, strcat('means_', properties{1}), 'A1');
    xlswrite(output_file, header, strcat('stds_', properties{1}), 'A1');
    xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{1}), 'A2');
    xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{1}), 'A2');
      subplot(2,1,2);

       [header, means_t, vars_t, stds_t] = single_csv_data(files{1}, num(1,3), num(1,4),properties{2},1, columns);
    cols = numel(header);
    means= zeros(files_num, cols);
    vars = zeros(files_num, cols);
    stds = zeros(files_num, cols);
    means(1,1:end) = means_t;
    vars(1,1:end) = vars_t;
    stds(1,1:end) = stds_t;
    for d = 2:files_num
        [header, means(d,1:end), vars(d,1:end), stds(d,1:end)] = single_csv_data(files{d}, num(d,3), num(d,4),properties{2},1,columns);
    end
     means(1:end,1)= num(1:end,1);
     vars(1:end,1)= num(1:end,1);
     stds(1:end,1)= num(1:end,1);


      %p = plot( fliplr(-(cols-4):-1).',(means(1:end,3:cols - 2)./repmat(means(1:end,cols-2),1,cols - 4)).', '.');
      %p = plot( (1:cols-3).',(means(1:end,3:cols - 1)).', '.');
      %rotate(p, [0 0 1], 270, [0 0 0]);
      %axis([1,10,0,5]);
      hold on;
      %for i = 1:files_num
          %plot((means(1:2:end, 3:cols - 1)).', 1:cols - 3, '.');
          %plot((means(2:2:end, 3:cols - 1)).', 1:cols - 3, 's','MarkerSize',3);
          layers = (4-cols:1:0);
          plot((means(1:end, 3:cols - 1)).', layers, '.');
          
          %e = errorbar((1:cols - 3).',(means(i,3:cols - 1)).',stds(i, 3:cols-1).', '.');
          %e.LineStyle = 'none';
          %set(e, 'MarkerFaceColor',colors(i, 1:3));
      %end
      title('Velocity Profile');
      xlabel(strcat('mean layer velocity'));
      ylabel('depth layer');
      
       furb_vel = 'u*log((0.5-x)/10)';
       
%        sta_v_fixed_mu_beta =  ...
%             strcat(2*rho*g* (num2str(sind(theta)) - st2_f.m_d*(cosd(theta)))/ ...
%                             (st2_f.b * 3*d*sqrt(rho_p/rho*g*cosd(theta))), ...
%                             '*(',num2str(depth^1.5),' - ((0.5-x)*',num2str(depth/10),').^1.5)' );
%        
%                      
%       
%       for theta_i = 1:knick
%         
%           i_fit = fit (layers.',means(theta_i, 3:cols - 1).', furb_vel);
%           i_plot = plot(i_fit,(means(1:end, 3:cols - 1)).', layers.');
%           set(i_plot(2), 'color', colors(theta_i,1:3));
%       end
      
      
      for theta_i=knick:numel(means(1:end,1))
          theta = means(theta_i,1);
          standa_velocity = strcat(num2str(2*rho*g),'* (',num2str(sind(theta)),' - m_d*', num2str(cosd(theta)), ...
                            ')/(b* ', num2str(3*d*sqrt(rho_p/rho*g*cosd(theta))),') * (', ...
                            num2str(depth^1.5),' - ((0.5-x)*',num2str(depth/10),').^1.5) ');
          i_fit = fit (layers.',means(theta_i, 3:cols - 1).', standa_velocity);
          %i_fit = fit (layers.',means(theta_i, 3:cols - 1).', sta_v_fixed_mu_beta);
          i_plot = plot(i_fit,(means(1:end, 3:cols - 1)).', layers.');
          set(i_plot(2), 'color', colors(theta_i+1,1:3));
      end
      lh = legend(arrayfun(@num2str, num(1:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
      lh = legend(arrayfun(@num2str, num(1:2:end,1), 'UniformOutput', false), 'location', 'southeast', 'FontSize',7);
         plot((means(1:end, 3:cols - 1)).', layers, '.');
      hold off;     
      saveas(gcf, strcat('../matlab_2/output/plot_',sub_name, '_knick_', num2str(knick)),'fig');
    %end

    xlswrite(output_file, header, strcat('means_', properties{2}), 'A1');
    xlswrite(output_file, header, strcat('stds_', properties{2}), 'A1');
    xlswrite(output_file, means(1:end,1:end), strcat('means_', properties{2}), 'A2');
    xlswrite(output_file, stds(1:end,1:end), strcat('stds_', properties{2}), 'A2');
    xlswrite(output_file, {'' 'K' 'Sc'  'm_s' 'm_d' 'beta' 'R rsquare' 'ST rsquare'}, 'fits', 'A1');
    xlswrite(output_file, {'' num2str(roer1_f.a) num2str(roer1_f.b)  num2str(st2_f.m_d) num2str(st2_f.b) num2str(roer1_goodness.rsquare) num2str(fern2_goodness.rsquare) num2str(st2_goodness.rsquare)}, 'fits', 'A2');
   
