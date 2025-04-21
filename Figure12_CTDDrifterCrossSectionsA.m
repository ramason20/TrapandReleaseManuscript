%Plots CTD Cast and Drifter Data
%Intended for Figure 12 Part A
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 6/28/22 
%Edited: 4/21/25
%------------------------------------------------------------------------

%% - Set Flags

isp = 0;    %Save plots? Y=1 N=0
dirname = 'F:\RAM_Research\Drifter_Data\2022June\CTD_Data\CrossSection_Figures\MABPOM_Figures\';

%% - Load In Data

load('F:\RAM_Research\Drifter_Data\2022June\CTD_Data\CTDData_2022Jun14.mat');

%% - Calculate Distance Between Casts

tran3_dist = zeros((length(tran3_lat)),1);

%Distance for Transect 3
for ilatlon = 1:(length(tran3_lat)-1)
    lat1 = tran3_lat(ilatlon)*(pi/180); lon1 = tran3_lon(ilatlon)*(pi/180);  
    lat2 = tran3_lat(ilatlon+1)*(pi/180); lon2 = tran3_lon(ilatlon+1)*(pi/180);
    distrad = 2*asin(sqrt((sin((lat1-lat2)/2))^2+cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2));
    tran3_dist(ilatlon+1) = distrad*(180/pi)*60*1.852; %In kilometers!
end

tran3_cumdist = cumsum(tran3_dist);

%% - Interplotate to an even depth grid for more even plotting!!

% Transect 3
Dq = 0:.05:21;
Vq3_Salt=zeros(length(Dq),size(tran3_Salt,2));
Vq3_dep=zeros(length(Dq),size(tran3_Salt,2));
Vq3_dis=zeros(length(Dq),size(tran3_Salt,2));
for iwd = 1:size(tran3_Salt,2)
    castd = tran3_Dep(:,iwd);   casts = tran3_Salt(:,iwd);
    castd(isnan(castd))=[];     casts(isnan(casts))=[];
    castd=cat(1,0,castd);       casts=cat(1,NaN,casts);
    
    Vq = interp1(castd,casts,Dq);
    Vq3_Salt(:,iwd) = Vq;
    Vq3_dep(:,iwd)=Dq;
    Vq3_dis(:,iwd)=tran3_cumdist(iwd);
end

%% - Pull Drifter Data and Find Nearest Point on the Cross-Section

load('F:\RAM_Research\Drifter_Data\2022June\Drifter_2022June14-15_Interp.mat');
itdrift3 = [59, 55, 51, 65, 56, 55, 62, 58, 61, 62];

t3dplat = [latint_2693514(itdrift3(1)),latint_2693519(itdrift3(2)),latint_4337850(itdrift3(3)),latint_4350132(itdrift3(4)),latint_4350393(itdrift3(5)),...
    latint_4351012(itdrift3(6)),latint_4351043(itdrift3(7)),latint_4351216(itdrift3(8)),latint_4354955(itdrift3(9)),latint_4355241(itdrift3(10))];
t3dplon = [lonint_2693514(itdrift3(1)),lonint_2693519(itdrift3(2)),lonint_4337850(itdrift3(3)),lonint_4350132(itdrift3(4)),lonint_4350393(itdrift3(5)),...
    lonint_4351012(itdrift3(6)),lonint_4351043(itdrift3(7)),lonint_4351216(itdrift3(8)),lonint_4354955(itdrift3(9)),lonint_4355241(itdrift3(10))];


%Finds the closest grid point (ix,iy)----------------------------------
satstepsize = 50;

    xsdrift = -t3dplon;    ysdrift = t3dplat;
    xr=linspace(tran3_lon(1),tran3_lon(2),satstepsize);
    yr=linspace(tran3_lat(1),tran3_lat(2),satstepsize);
    for idrift = 2:(length(xsdrift)-1)
        xfp=linspace(tran3_lon(idrift),tran3_lon(idrift+1),satstepsize);
        yfp=linspace(tran3_lat(idrift),tran3_lat(idrift+1),satstepsize);
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
        lat1 = tran3_lat(1)*(pi/180);   lon1 = tran3_lon(1)*(pi/180);  
        lat2 = iys(ip2)*(pi/180);       lon2 = ixs(ip2)*(pi/180);
        distrad = 2*asin(sqrt((sin((lat1-lat2)/2))^2+cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2));
        drift3_dist(ip2) = distrad*(180/pi)*60*1.852; %In kilometers!
    end

%% - Plot Data

saltmax = 28;
saltmin = 22;

%Transect 3
cmap = cmocean('haline');
figure
set(gcf,'Position',[150   150   1400   400])
%pcolor(tran3_cumdist,-tran3_Dep,tran3_Salt);
pcolor(Vq3_dis,-Vq3_dep,Vq3_Salt);
shading('flat');
hold on
[Conman,hlab] = contour(Vq3_dis,-Vq3_dep,Vq3_Salt,[saltmin:.25:saltmax],'k-');
plot(drift3_dist,zeros(1,10),'rd','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
%vlab = saltmin:0.5:saltmax;
vlab = saltmin:saltmax;
clabel(Conman,hlab,vlab)
grid on
c1 = colorbar;
colormap(cmap)
caxis([saltmin saltmax])
c1.Label.String = 'Salinity (PSU)';
xlabel('Cross-Channel Distance (km)');
xlim([0 11])
ylabel('Depth (m)')
title('Lower Bay Spring Tide')
ylim([-21 0.5])
set(gca,'FontSize',16)
print([dirname,'Tran3_2022Jun14'],'-dpng');
