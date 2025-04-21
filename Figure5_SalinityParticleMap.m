%%%Creates Salinity Plots Based Off the Delaware Bay ROMS Salinty Field
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 10/18/24 
%Edited: 4/21/25
%------------------------------------------------------------------------

%% - Flags

tic
springtide = 1;
isv = 0;
tidavg = 1;

%% - Load Data

if springtide == 1
    load('\SpringParticles\PartTrack\11.mat');
    partx = xs(84,:);    party = ys(84,:);      %High Slack
    %partx = xs(120,:);   party = ys(120,:);     %Max Flood
    %partx = xs(60,:);    party = ys(60,:);      %Max Ebb
    %partx = xs(36,:);    party = ys(36,:);      %Low Slack
    saltmonth = '\OceanModels\ocean_his_0010.nc';
    if tidavg == 1
        salttime = 509:558;
    else
        salttime = 573; %HS
        %salttime = 576; %Max Flood
        %salttime = 570; %Max Ebb
    end
else
    load('\NeapParticles\PartTrack\11.mat');
    partx = xs(90,:);    party = ys(90,:);      %High Slack
    %partx = xs(120,:);   party = ys(120,:);     %Max Flood
    %partx = xs(60,:);    party = ys(60,:);      %Max Ebb
    %partx = xs(36,:);    party = ys(36,:);      %Low Slack
    saltmonth = '\OceanModels\ocean_his_0011.nc';
    if tidavg == 1
        salttime = 3:52;
    else
        salttime = 29;
    end
end

saltdata = nc_varget(saltmonth,'salt');     %Load salinity data
saltdata = squeeze(saltdata(salttime,10,:,:));     %Squeeze salinity data down to the surface layer
if tidavg == 1
    saltdata = mean(saltdata,1);
    saltdata = squeeze(saltdata(1,:,:));
end

lat_rho = nc_varget(saltmonth,'lat_rho');
lon_rho = nc_varget(saltmonth,'lon_rho');

ot=nc_varget(saltmonth,'ocean_time');       %Load the date/time parameter from the netcdf file
t0=datenum(2016,9,1,0,0,0);                 %Initial time for the ROMS run (for 2008 - 2008,1,1,0,0,0)
otc = (ot/3600/24)+t0;                      %Converts the ot time steps to calendar day/hour/minute

load('\ParticleTracks\StartUpFile_Step50m.mat')  %Load 
clear excelpart     %Unneeded large file for plotting salinity (needed for particle initialization in particle tracking)

%UD Cross
x1 = 216.753;   y1 = 319.825;   x2 = 208.949;   y2 = 331.991;
alongslope = (y2-y1)/(x2-x1);   crossslope = -(1/alongslope);
brandyx = 216.011;              brandyy = 327.993;
uppointx = 221.5;               uppointy = (crossslope*(uppointx-brandyx))+brandyy;
dnpointx = 208.5;               dnpointy = (crossslope*(dnpointx-brandyx))+brandyy;
crossx = [dnpointx uppointx];   crossy = [dnpointy uppointy];

%Creation of Kill Line
xDE = 220;  yDE = 298;  %First values used: x=215 y=295
xNJ = 250;  yNJ = 334;  %First values used: x=245 y=330
numopoints = 1000;
outxline=linspace(xDE,xNJ,numopoints);
outyline=linspace(yDE,yNJ,numopoints);

%Convert from xr/yr into lon/lat
xrs = reshape(xr,[300*150, 1]);
yrs = reshape(yr,[300*150, 1]);
lonrs = reshape(lon_rho,[300*150, 1]);
latrs = reshape(lat_rho,[300*150, 1]);
Fx = scatteredInterpolant(xrs,yrs,lonrs);
Fy = scatteredInterpolant(xrs,yrs,latrs);

xpart = Fx(partx,partx);        ypart = Fy(partx,party);
xcross = Fx(crossx,crossy);     ycross = Fy(crossx,crossy);
xout = Fx(outxline,outyline);   yout = Fy(outxline,outyline);

%% - Plot Cross-Section

fig1=figure;                                    % Used for plotting main figure
set(gcf,'PaperSize',[13.33/2 7.5])              % Sets scale of figure so it does not become distorted
set(gcf,'Color','w')                            % Sets plot background color to be black
set(gcf, 'InvertHardcopy', 'off')               % Tells Matlab to save the figure background color and not to change it to white
set(gcf,'Position',[378   302   550   653])     % Sets location and size of the plot
cnt = 1000;                                     % Imagine number counting for saving figures
cmap = cmocean('ice');    

%Plots Bathymetry Data
bathy(250:300,1:60)=nan;                     %Eliminate the artifical channel data that appears in the ploting area
pcolor(lon_rho,lat_rho,-bathy); shading('flat')
hold on
plot(xcross,ycross,'r-')
plot(xout,yout,'r--')
c1 = colorbar;
colormap(cmap)
caxis([-30 0])
c1.Label.String = 'Depth (m)';         %Colorbar label
set(gca,'DataAspectRatio',[1 1 1])          %Adjusts aspect ration of figure
set(gca,'XLim',[-75.5  -74.8],'YLim', [38.4  39.45])
set(gca,'FontSize',12)                      %Adjusts figure font size to be readable
xlabel('Longitude')
ylabel('Latitude')
grid on

if isv == 1
    savdir = '\savloc\';
    print([savdir,'BathyFigure'],'-dpng');
end

%% - Plot the Data

cmapX = cool(60); cmapY = spring(30); cmapZ = flipud(parula(15));
cmap = [cmapX;cmapY;cmapZ];

fig1=figure;                                    % Used for plotting main figure
set(gcf,'PaperSize',[13.33/2 7.5])              % Sets scale of figure so it does not become distorted
set(gcf,'Color','w')                            % Sets plot background color to be black
set(gcf, 'InvertHardcopy', 'off')               % Tells Matlab to save the figure background color and not to change it to white
set(gcf,'Position',[378   302   550   653])     % Sets location and size of the plot
cnt = 1000;                                     % Imagine number counting for saving figures
    
%Plots Bathymetry Data
salthour = saltdata;                            %Squeeze surface salinity down to current time step
salthour(250:300,1:60)=nan;                     %Eliminate the artifical channel data that appears in the ploting area
pcolor(lon_rho,lat_rho,salthour); shading('flat')
hold on
plot(xcross,ycross,'r-')
plot(xpart,ypart,'k.','MarkerSize',0.05)
c1 = colorbar;
colormap(cmap)
caxis([0 35])
c1.Label.String = 'Salinity (PSU)';         %Colorbar label
set(gca,'DataAspectRatio',[1 1 1])          %Adjusts aspect ration of figure
%set(gca,'XColor',0.7*[1 1 1],'YColor',0.7*[1 1 1])  %Changes color of x/y axis text/labels
set(gca,'XLim',[-75.5  -74.8],'YLim', [38.4  39.45])
set(gca,'FontSize',12)                      %Adjusts figure font size to be readable
xlabel('x (km)')    
xlabel('Longitude')
ylabel('Latitude')
grid on

if isv == 1
    savdir = '\savloc\';
    if tidavg == 1
        print([savdir,'PaperFigure_TidalAvgSalt_wcross'],'-dpng');
        %print([savdir,'PaperFigure_TidalAvgSalt_MaxFloodPart_wcross'],'-dpng');
        %print([savdir,'PaperFigure_TidalAvgSalt_MaxEbbPart_wcross'],'-dpng');
        %print([savdir,'PaperFigure_TidalAvgSalt_LowSlackPart_wcross'],'-dpng');
    else
        print([savdir,'PaperFigure_Hour127_wcross'],'-dpng');
    end
end

toc

