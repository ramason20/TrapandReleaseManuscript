%Plots CTD Cast and Drifter Data
%Intended for Figure 12 Parts B and C
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 6/28/24 
%Edited: 4/21/25
%------------------------------------------------------------------------

%% - Set-Up Flags
spring = 1;
neap = 0;
isv = 1;
dirname = 'C:\Savloc\';

%% - Load CTD Data
if spring == 1
    load('C:\UpperBaySpringTideCTD.mat');
    tranno = 2;
    transectfull = 1:8;
elseif neap == 1
    load('C:\UpperBayNeapTideCTD.mat');
    tranno = 16;
    transectfull = 1:54;
end

%Find Distance Between CTD Casts
    tottran = max(survnum);
    tranf = transectfull(1);    trane = transectfull(end);
    TranSalt = squeeze(salt(:,tranf:trane));
        
    lattran = squeeze(lat(tranf:trane));    lontran = squeeze(lon(tranf:trane));
    TranDist = zeros(1,length(lattran));
    for ilat = 2:length(lattran)
        lat1 = lattran(1)*(pi/180);      lon1 = lontran(1)*(pi/180);  
        lat2 = lattran(ilat)*(pi/180);   lon2 = lontran(ilat)*(pi/180);
        distrad = 2*asin(sqrt((sin((lat1-lat2)/2))^2+cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2));
        TranDist((ilat)) = distrad*(180/pi)*60*1.852*1000;
    end
    TranDistFull = ones(size(TranSalt,1),size(TranSalt,2)).*TranDist;
        
    depthtran = squeeze(depth(tranf:trane));
    TranDepth = ones(size(TranSalt,1),size(TranSalt,2));
    depthar = linspace(0,max(depth),size(TranSalt,1));
    for id = 1:size(TranSalt,2)
        TranDepth(:,id) = (-1.*depthar)';
    end

%% - Load Drifter Data

if spring == 1
    load('C:\2023Aug01Release.mat');
    t0EDT = datetime(2023,08,01,10,00,0);
    tendEDT = datetime(2023,08,02,11,00,0);
    tdrift = t0EDT:minutes(5):tendEDT;
    castpoint = 57;
    drifterloncat = cat(2,drift1_lon,drift2_lon,drift3_lon,drift4_lon,drift5_lon,drift6_lon);
    drifterlatcat = cat(2,drift1_lat,drift2_lat,drift3_lat,drift4_lat,drift5_lat,drift6_lat);
    drifterlon = drifterloncat(castpoint,:);
    drifterlat = drifterlatcat(castpoint,:);
elseif neap == 1
    load('C:\2023Jul27Release.mat');
    t0EDT = datetime(2023,07,27,9,00,0);
    tendEDT = datetime(2023,07,31,11,00,0);
    tdrift = t0EDT:minutes(5):tendEDT;
    castpoint = 121;
    drifterloncat = cat(2,drift1_lon,drift2_lon,drift3_lon,drift5_lon,drift6_lon);
    drifterlatcat = cat(2,drift1_lat,drift2_lat,drift3_lat,drift5_lat,drift6_lat);
    drifterlon = drifterloncat(castpoint,:);
    drifterlat = drifterlatcat(castpoint,:);
end

%Finds the closest grid point (ix,iy)----------------------------------

    satstepsize = 25;
    xsdrift = drifterlon;    ysdrift = drifterlat;
    xr=linspace(lontran(1),lontran(2),satstepsize);
    yr=linspace(lattran(1),lattran(2),satstepsize);
    for idrift = 2:(length(xsdrift)-1)
        xfp=linspace(lontran(idrift),lontran(idrift+1),satstepsize);
        yfp=linspace(lattran(idrift),lattran(idrift+1),satstepsize);
        xr = cat(2,xr,xfp(2:end));
        yr = cat(2,yr,yfp(2:end));
    end

    np = length(xsdrift);
    ixs = zeros(1,np); iys=zeros(1,np);     % Grid Location Array
    for ip=1:np
        x1 = xr(1);     y1 = yr(1);
        x2 = xr(end);   y2 = yr(end);
        x3=xsdrift(1,ip); y3=ysdrift(1,ip);
        k = (((y2-y1)*(x3-x1))-((x2-x1)*(y3-y1)))/(((y2-y1)^2)+((x2-x1)^2));
        ixs(ip) = x3-(k*(y2-y1));
        iys(ip) = y3+(k*(x2-x1));
    end
    
    for ip2 = 1:np
        lat1 = lattran(1)*(pi/180);     lon1 = lontran(1)*(pi/180);  
        lat2 = iys(ip2)*(pi/180);       lon2 = ixs(ip2)*(pi/180);
        distrad = 2*asin(sqrt((sin((lat1-lat2)/2))^2+cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2));
        drift_dist(ip2) = distrad*(180/pi)*60*1.852; %In kilometers!
    end
    


%% - Plot data

        if neap == 1
            %cmap = cmocean('balance',40);
            cmap = cmocean('haline');
        else
            %cmap = cmocean('balance',56);
            cmap = cmocean('haline');
        end
       
        figure
        set(gcf,'Position',[150   150   1400   400])
        pcolor((TranDistFull./1000),TranDepth,TranSalt)
        hold on
        contour((TranDistFull./1000),TranDepth,TranSalt,[4:20],'k','ShowText','on')
        plot((TranDist./1000),-depthtran','k-')
        plot(drift_dist,zeros(1,np),'rd','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
        shading('flat')
        xlabel('Cross Channel Distance (km)')
        xlim([0 3])
        ylabel('Depth (m)')
        ylim([-20 0.5])
        if spring == 1
            title('Upper Bay Spring Tide')
        elseif neap == 1
            title('Neap Tide')
        end
        grid on
        c1 = colorbar;
        colormap(cmap)
        c1.Label.String = 'Salinity (PSU)';
        if spring == 1
            caxis([10 13])
        elseif neap == 1
            caxis([9 19])
        end
        set(gca,'FontSize',16)
        if isv == 1 
            if spring == 1
                print([dirname,'SpringTide'],'-dpng');
            elseif neap == 1
                print([dirname,'NeapTide'],'-dpng');
            end
        else
            pause(1)
        end
        %close all
