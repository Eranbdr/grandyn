function [p_fft, p_raw  ] = csv_fft_single_fig(file_name, dt)
 %[csv_dat, header] = xlsread(file_name, strcat(first_col, num2str(1),':',last_col, num2str(1)));    
%[csv_data] = xlsread(file_name, strcat(first_col, num2str(first_id),':',last_col, num2str(last_id-1)));
[csv_data] = csvread(file_name);
    set(0,'DefaultAxesColorOrder',[0.5 0.5 0.5; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 0 0 0; 1 1 0; 0.2 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2; 0.2 0.8 0.8; 0.8 0.2 0.2; 0.8 0.2 0.8; 0.8 0.8 0.2; 0.8 0.8 0.8]);
    colors = get(0, 'DefaultAxesColorOrder');
    
    data = csv_data(5000000:35000000);%csv_data(first_id:last_id-1,1:cols);
    [rows, cols] = size(data);
    L = length(data);
    TotalTime = L*dt;
    Fs = L/TotalTime;
    t = dt*(0:1:L-1);
    n = ceil(log2(length(data)));
    f=Fs*(0:(2^(n-1)-1))/2^n;
    period = 1./f;
    
    smooth_data = zeros(L, cols);
    for i = 1:cols
    
        smooth_data(1:end,i) = data(1:end,i) - mean(data(1:end,i));
        %smoothflux = filter(ones(1,5)/5,1,totalflux);
        %smoothflux=detrend(totalflux);
    end
        %totalflux= totalflux(1:100:end);


        

        %plot(t,smooth_data,'bx-');
        %NFFT = 2^nextpow2(L);

        Y = fft(smooth_data, 2^n);
        Pyy = Y.*conj(Y)/2^n;
        
        %converting frequncy to period.
    %figure('name', strrep(strrep(file_name,'_', ' '), '.csv', ''));
    
    title(strcat('fft ', strrep(file_name, '_', ' ')));
   %figure;
    subplot(2,1,1);
    %hold on;   
        %figure(2)
        p_fft = plot(period(2:end),Pyy(2:2^(n-1),1:end),'-');
        %figure(3)
        %plot(f(2:end),Pyy(2:2^(n-1)),'rx-');
    
    
    subplot(2,1,2);
    p_raw = plot((dt:dt:rows*dt), data(1:end,1:end));
    
    ylabel('[velocity / flux]');
    xlabel('time steps');
    %legend(strrep (header(1:cols), '_',' '), 'location', 'west');
    saveas(gcf, strcat('fft_', strrep(strrep(file_name, '_', ' '), '.csv', '.fig')));

end

