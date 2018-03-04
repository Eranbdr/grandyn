function RestartData = MyReadRestart(name, inum, fnum, interval, orig_height, layers)
i = inum;
if i < 10
    filename = strcat('restart000',num2str(i));
elseif i < 100
    filename = strcat('restart00',num2str(i));
elseif i < 1000
    filename = strcat('restart0',num2str(i));
else
    filename = strcat('restart',num2str(i));
end
fid = fopen(filename,'r');
garbage = fread(fid,11,'schar');
garbage = fread(fid,2,'int32');
n = fread(fid,1,'int32');
total_flux = 0;
velocities = zeros(layers,n);
RestartData = zeros((fnum-inum)/interval,n,12);
for i = inum:fnum
    if i < 10
        filename = strcat('restart000',num2str(i));
    elseif i < 100
        filename = strcat('restart00',num2str(i));
    elseif i < 1000
        filename = strcat('restart0',num2str(i));
    else
        filename = strcat('restart',num2str(i));
    end
    fid = fopen(filename,'r');
    garbage = fread(fid,11,'schar');
    garbage = fread(fid,2,'int32');
    n = fread(fid,1,'int32');
    garbage = fread(fid,4,'int32');
    garbage = fread(fid,1,'float64');
    ybop = fread(fid,1,'float64');
    ytop = fread(fid,1,'float64');
    xleft = fread(fid,1,'float64');
    xright = fread(fid,1,'float64');
    garbage = fread(fid,1,'float64');
    tau = fread(fid,1,'float64');
    tauadd =  fread(fid,1,'float64');
    garbage = fread(fid,1,'float64');
    DataGrain = zeros(n,12);
    for j = 1:n
        rx = fread(fid,1,'float64');
        ry = fread(fid,1,'float64');
        radius = fread(fid,1,'float64');
        rnt = fread(fid,1,'float64');
        vx = fread(fid,1,'float64');
        vy = fread(fid,1,'float64');
        w = fread(fid,1,'float64');
        fx = fread(fid,1,'float64');
        fy = fread(fid,1,'float64');
        tq = fread(fid,1,'float64');
        bdgrain = fread(fid,1,'int32');
        gtype = fread(fid,1,'int32');
        DataGrain(j,:) = [rx ry radius rnt vx vy w fx fy tq bdgrain gtype];
    end
    fclose(fid);
    RestartData(i-inum+1,:,:) = DataGrain(:,:);
    
end
    
    
        