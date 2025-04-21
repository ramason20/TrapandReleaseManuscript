%%% Plot Code of Cross-Section for the Idealized Model
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 9/19/24 
%Edited: 4/21/25
%------------------------------------------------------------------------

%%%Plot Coriolis Idealized Model Data

%Non-Rotational Trapping
load('\data\rk_suvw_nx151B10000h12Km275g19097sx200it2000000.mat');
savdir = '\savloc\';

%Non-Rotational Flushing
% load('\data\rk_suvw_nx151B10000h12Km50g099994sx200it2000000.mat');
% savdir = 'savloc\';

%Rotational Trapping
% load('\data\rk_suvw_nx151B10000h12Km275g19097sx200it20000000.mat');
% savdir = '\savloc\';

%Rotational Flushing
% load('\data\rk_suvw_nx151B10000h12Km50g099994sx200it20000000.mat');
% savdir = '\savloc\';

%% - Put v and w on the rho grid

varray = zeros(size(u,1),size(u,2));
for iv = 1:(size(v,1)-1)
    for i = 1:size(v,2)
        va1 = v(iv,i);     va2 = v(iv+1,i);
        varray(iv,i)=(va1+va2)/2;
    end
end

warray = zeros(size(u,1),size(u,2));
for iw = 1:size(w,1)
    for ii = 1:(size(w,2)-1)
        wa1 = w(iw,ii);     wa2 = w(iw,ii+1);
        warray(iw,ii)=(wa1+wa2)/2;
    end
end

%% - Set-Up Plot Variables

depthplotr = ones(size(u,1),size(u,2)).*((-(nz-1):0).*dz);
widthplotr = ones(size(u,1),size(u,2)).*(((0:(ny-1)).*dy)');
depthplotv = ones(size(v,1),size(v,2)).*((-(nz-1):0).*dz);
widthplotv = ones(size(v,1),size(v,2)).*(((0:(ny)).*dy)');
widthplotrn = (widthplotr-5000)./1000;
widthplotvn = (widthplotv-5000)./1000;

%% - Plot Across-Contours over Along-Channel Vel

fig=figure;
set(gcf,'Position',[150   50   1400   500])
pcolor(widthplotrn',depthplotr',-u'); shading('flat')
xlim([-5 5])
xlabel('Distance (km)','FontSize',20)
ylabel('Depth (m)','FontSize',20)
ax = gca; 
ax.FontSize = 20;
c1 = colorbar;
colormap(cmap)
%caxis([-0.06 0.06]);
%caxis([-0.125 0.125]); %Flushing
caxis([-0.025 0.025]);  %Trapping
%caxis([-0.9 0.9]);
c1.Label.String = 'Along-Channel Velocity (m s^{-1})';
hold on
%conval = -1:0.0025:1;
conval = -1:0.01:1;
title('Along-Channel Velocity with Across-Channel Velocity Contours')
[C,h]=contour(widthplotrn',depthplotr',varray',conval,'k-','ShowText','on');
clabel(C,h,'FontSize',16,'Color','black')
grid on
print([savdir,'AlongAcrossPlot'],'-dpng');

%% - Plot Across-Contours over Salinity

figure
set(gcf,'Position',[150   50   1400   500])
pcolor(widthplotrn',depthplotr',s'); shading('flat')
xlim([-5 5])
xlabel('Distance (km)','FontSize',20)
ylabel('Depth (m)','FontSize',20)
ax = gca; 
ax.FontSize = 20;
c1 = colorbar;
colormap(cmap)
caxis([-2 2]);  %Non-Rotating
%caxis([-3 3]);  %Rotational
c1.Label.String = 'Salinity (PSU)';
hold on
%conval = -1:0.0025:1;
conval = -1:0.01:1;
title('Salinity with Across-Channel Velocity Contours')
[C,h] = contour(widthplotrn',depthplotr',varray',conval,'k-','ShowText','on');
clabel(C,h,'FontSize',16,'Color','black')
grid on
print([savdir,'SaltAcrossPlot'],'-dpng');
