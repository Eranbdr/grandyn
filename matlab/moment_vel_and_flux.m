function[] = moment_vel_and_flux(dir, file_name, num_of_layers, roer_step, bag_step, start, end_, spread)
%moment_vel('../matlab/', 'velocity_profile_random_M1_15_19gy_19gx_12000WL_5_restart', 
%                      20,1417300000/50000,803400000/50000)

% start = 1000;
% end_ = 500;
% spread = 500;


[flux_data, header] = xlsread(strcat(dir, 'flux_', file_name,'.csv'));
[vel_data, header] = xlsread(strcat(dir, 'velocity_profile_', file_name,'.csv'));
 fig_H  = figure('name', file_name);
         
       %legend off;
          subplot(2,1,1);
        hold on;  
        title('Simuation at Transition Angle','FontSize',14); 
    xlabel('Time Step','FontSize',12); 
    ylabel('Flux (non-dimensional)','FontSize',12); 
    %legend('location', 'northwest');
    
    plot(flux_data(:,1), flux_data(:,num_of_layers + 3),'b', 'LineWidth', 2);
    roer_flux_min = 0;
    bag_flux_max = 0;
    
    if roer_step < 10
     

      % sort according to total flux
       [max_flux, max_id]= max(flux_data(start:end-end_,num_of_layers+3));
       [min_flux, min_id]= min(flux_data(start:end-end_,num_of_layers+3));
      
      %circles of max and min flux values
       scatter(flux_data(start+min_id,1), min_flux,'r','LineWidth',3);
      scatter(flux_data(start+max_id,1), max_flux,'k', 'LineWidth',3);
      if roer_step > 1
      % sort according to total flux
      
       [max_flux1, max_id1]= max(flux_data(start:max_id+start-spread,num_of_layers+3));
       [min_flux1, min_id1]= min(flux_data(start:min_id+start-spread,num_of_layers+3));
      
      %circles of max and min flux values
       scatter(flux_data(start+min_id1,1), min_flux1,'m','LineWidth',3);
      scatter(flux_data(start+max_id1,1), max_flux1,'g', 'LineWidth',3);
      if roer_step > 2
      % sort according to total flux
       [min_flux2, min_id2]= min(flux_data(min_id+start+spread:end-end_,num_of_layers+3));
       [max_flux2, max_id2]= max(flux_data(max_id+start+spread:end-end_,num_of_layers+3));
      
      %circles of max and min flux values
       scatter(flux_data(min_id+start+spread+min_id2,1), min_flux2,'c','LineWidth',3);
      scatter(flux_data(max_id+start+spread+max_id2,1), max_flux2,'y', 'LineWidth',3);
      end
      end
    else
        scatter(flux_data(roer_step,1), flux_data(roer_step,num_of_layers + 3),'r')
        scatter(flux_data(bag_step,1), flux_data(bag_step,num_of_layers + 3),'k')
    
    end
    
    
    
    sub_ax = subplot(2,1,2);
  title ('Instantaneous Velocity Profile','FontSize',14); 
  h=68.8;
  depth =  [ (h  -h/(2*num_of_layers): -h/(num_of_layers) : h/(2*num_of_layers))].';
  %view(sub_ax,[90 -90]);
        hold on;
        
        if roer_step < 10
            
            % how many profiles to draw?
          plot(depth, vel_data(min_id+start,3:num_of_layers + 2),'.r');
          plot(depth, vel_data(max_id+start,3:num_of_layers + 2),'.k');
  if roer_step > 1
          if  numel(min_id1) == 1; plot(depth, vel_data(min_id1+start,3:num_of_layers + 2),'.m');end;
           if  numel(max_id1) == 1;plot(depth, vel_data(max_id1+start,3:num_of_layers + 2),'.g');end;
         if roer_step > 2 
           if  numel(min_id2) == 1;plot(depth, vel_data(max_id+start+spread+min_id2,3:num_of_layers + 2),'.c');end;
           if  numel(max_id2) == 1;plot(depth, vel_data(max_id+start+spread+max_id2,3:num_of_layers + 2),'.y');end;
         end
  end
        else
            
            % single entry
            plot(depth, vel_data(roer_step,3:num_of_layers + 2),'.r')
             plot(depth, vel_data(bag_step,3:num_of_layers + 2),'.k')
        end
        
        
    view(sub_ax,[90 90]);
      ylabel('Velocity (non-dimensional)','FontSize',12); 
    xlabel('Depth (# grains)','FontSize',12); 
    saveas(gcf, strcat('output/instant/inst_',strrep(file_name, '.','-')));
end
    