% data reading
%ERA5
era5path='G:\ERA5daily\';
monvar={'01','02','03','04','05','06','07','08','09','10','11','12'};
rhobs36y=[];tdobs36y=[];
for y=1:36
    year=y+1978;
    disp(year)
    tdy=[];d2y=[];   
    for mon=1:12
%         ncinfo([era5path,'download_daily_mean_2m_dewpoint_temperature_',num2str(year),'_0',num2str(mon),'.nc']);
        tdm=single(ncread([era5path,'download_daily_mean_2m_temperature_',num2str(year),'_',monvar{mon},'.nc'],'t2m')); 
%         psm=single(ncread([era5path,'download_daily_mean_surface_pressure_',num2str(year),'_',monvar{mon},'.nc'],'sp'));  
        d2m=single(ncread([era5path,'download_daily_mean_2m_dewpoint_temperature_',num2str(year),'_',monvar{mon},'.nc'],'d2m')); 
        tdy=cat(3,tdy,tdm);   
%         psy=cat(3,psy,psm);    
        d2y=cat(3,d2y,d2m); 
    end
    clear tdm psm d2m mon
    d2y=d2y-273.15; tdy=tdy-273.15; 
    e2my=6.11*10.^((7.63*d2y)./(241.9+d2y));
    es2my=6.11*10.^((7.63*tdy)./(241.9+tdy)); 
    rhy=e2my./es2my*100; 
    rhobs36y=cat(3,rhobs36y,rhy);    
    tdobs36y=cat(3,tdobs36y,tdy);
    twy=tdy.*atan(0.151977*(rhy+8.313659).^0.5)+atan(tdy+rhy)-atan(rhy-1.676331)+0.00391838*(rhy.^0.5).*atan(0.023101*rhy)-4.686035;
    thy=twy+4.5*(1-(rhy/100).^2); 
    thy=rot90(thy);twy=rot90(twy);    
    eval(strcat("save",{32},"G:\global\tw",int2str(year),".mat",{32},"twy"));
    eval(strcat("save",{32},"G:\global\th",int2str(year),".mat",{32},"thy"));
end
save G:\global\rhobs36y.mat rhobs36y
save G:\global\tdobs36y.mat tdobs36y

% pr
path='G:\cpc\cpc global precipitation£¨79-20nc£©\mat\';
for y=1:36
    year=y+1978;
    disp(year)
    pry0=single(importdata([path,'cpc',int2str(year),'.mat']));
    pry0(pry0<0)=nan;  
    clear pry
    for lat=1:180 
        for lon=1:360 
            a=pry0(lat*2-1:lat*2,lon*2-1:lon*2,:);
            pry(lat,lon,:)=nanmean(a,[1,2]);
        end        
    end   
    eval(strcat("save",{32},"G:\global\pr",int2str(year),".mat",{32},"pry"));
end

%CMIP6
gcms={'CNRM-CM6-1','EC-Earth3-Veg','KACE-1-0-G','MPI-ESM1-2-HR','MRI-ESM2-0','NorESM2-MM'};
ssp={'ssp126';'ssp245';'ssp370';'ssp585'};
for n=1:6
    disp(gcms{n})
    pathmod=['H:\',gcms{n},'\mat\'];    
    new_folder = ['G:/',gcms{n},'/mat/th/'];    mkdir(new_folder);  
    new_folder = ['G:/',gcms{n},'/mat/tw/'];    mkdir(new_folder); 
    clear new_folder
    filestdhis=dir([pathmod,'tas_*historical*.mat']);
    filesrhhis=dir([pathmod,'hurs_*historical*.mat']);
    for y=1:36       
        tdy=rot90(single(importdata([pathmod,filestdhis(1).name(1:end-8),int2str(y+1978),'.mat'])))-273.15; 
        rhy=rot90(single(importdata([pathmod,filesrhhis(1).name(1:end-8),int2str(y+1978),'.mat']))); 
        twy=tdy.*atan(0.151977*(rhy+8.313659).^0.5)+atan(tdy+rhy)-atan(rhy-1.676331)+0.00391838*(rhy.^0.5).*atan(0.023101*rhy)-4.686035; 
        thy=twy+4.5*(1-(rhy/100).^2);   
        a=strcat('save',{32},'G:\',gcms{n},'\mat\tw\','tw',int2str(y+1978),'.mat',{32},'twy');
        eval(a{1});       
        b=strcat('save',{32},'G:\',gcms{n},'\mat\th\','th',int2str(y+1978),'.mat',{32},'thy');
        eval(b{1});      
        clear a b twy thy tdy rhy
    end
    for m=1:4
        disp(ssp{m})
        filestdssp=dir([pathmod,'tas_*',ssp{m},'*.mat']);
        filesrhssp=dir([pathmod,'hurs_*',ssp{m},'*.mat']);
        for y=1:86
            tdy=rot90(single(importdata([pathmod,filestdssp(1).name(1:end-8),int2str(y+2014),'.mat'])))-273.15; 
            rhy=rot90(single(importdata([pathmod,filesrhssp(1).name(1:end-8),int2str(y+2014),'.mat'])));
            twy=tdy.*atan(0.151977*(rhy+8.313659).^0.5)+atan(tdy+rhy)-atan(rhy-1.676331)+0.00391838*(rhy.^0.5).*atan(0.023101*rhy)-4.686035; 
            thy=twy+4.5*(1-(rhy/100).^2); 
            a=strcat('save',{32},'G:\',gcms{n},'\mat\tw\','tw',ssp{m},'_',int2str(y+2014),".mat",{32},'twy'); 
            b=strcat('save',{32},'G:\',gcms{n},'\mat\th\','th',ssp{m},'_',int2str(y+2014),".mat",{32},'thy');
            eval(a{1}); eval(b{1}); 
            clear a b twy thy tdy rhy
        end
    end   
end

% gpp
gcms={'CNRM-CM6-1','EC-Earth3-Veg','KACE-1-0-G','MPI-ESM1-2-HR','MRI-ESM2-0','NorESM2-MM'};
ssp={'ssp126';'ssp245';'ssp370';'ssp585'};
P_lonArray=single(0.5:1:359.5); P_latArray=single(89.5:-1:-89.5); [lon1,lat1]=meshgrid(P_lonArray,P_latArray);
clear P_lonArray P_latArray
for n=1:6
    ncfiles=dir(['H:\',gcms{n},'\gpp\gpp*historical*.nc']);
    lonnc1=single(ncread(strcat('H:\',gcms{n},'\gpp\',ncfiles(1).name),'lon'));
    latnc1=single(flipud(ncread(strcat('H:\',gcms{n},'\gpp\',ncfiles(1).name),'lat')));
    [lonnc,latnc]=meshgrid(lonnc1,latnc1);
    clear lonnc1 latnc1      
    if ~isempty(ncfiles)
        gppy_1_36y=[];
        for i=1:length(ncfiles)
            gppy1=rot90(ncread(['H:\',gcms{n},'\gpp\',ncfiles(i).name],'gpp'));
            gppy=change0_360to_180_180(gppy1);
            clear gppy1 
            for day=1:size(gppy,3)
                gppy_1(:,:,day)=(interp2(lonnc,latnc,gppy(:,:,day),lon1,lat1,'linear',0)).*BD;
            end
            gppy_1_36y=cat(3,gppy_1_36y,gppy_1);
            clear gppy_1
        end
        if n==6            
            gppy_1_36y(:,:,1:9*12)=[];
        elseif n==5
            gppy_1_36y(:,:,1:129*12)=[];             
        elseif n==4
            gppy_1_36y(:,:,1:4*12)=[];
        end
        eval(strcat('save',{32} ,"'H:\",gcms{n},"\gpp\gppy_1_36y.mat'",{32},'gppy_1_36y')); 
        clear ncfiles i
    end
    for m=1:4        
        ncfiles=dir(['H:\',gcms{n},'\gpp\gpp*',ssp{m},'*.nc']);
        if isempty(ncfiles)
            continue;
        else
            gppy_1_86y=[];
            for i=1:length(ncfiles)
                gppy1=rot90(ncread(['H:\',gcms{n},'\gpp\',ncfiles(i).name],'gpp'));
                gppy=change0_360to_180_180(gppy1);
                clear gppy1 
                for day=1:size(gppy,3)
                    gppy_1(:,:,day)=(interp2(lonnc,latnc,gppy(:,:,day),lon1,lat1,'linear',0)).*BD;
                end
                gppy_1_86y=cat(3,gppy_1_86y,gppy_1);
                clear gppy_1
            end
            eval(strcat('save',{32} ,"'H:\",gcms{n},'\gpp\gppy_1_86y',ssp{m},".mat'",{32},'gppy_1_86y'));
            clear ncfiles i
        end
    end 
end
