%%% UDel Model Rotate into the DBOFS Cross-Section Projection
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 9/19/24 
%Edited: 4/21/25
%------------------------------------------------------------------------

% - Last Edit = Initial Creation
tic

%% - Flags

spring = 1;
salt=1;
anglecheck = 0;
isv = 1;
dirname = '\dirname\';

%% - Load in UD Model Data

load('C\dirname\StartUpFile.mat');

ix = 204;   iy = 50;
cmap = cmocean('ice');

figure
pcolor(xr,yr,bathy)
shading('flat')
hold on
plot(xr(ix,iy),yr(ix,iy),'rp')
plot([216.753 208.949],[319.825 331.991],'r-')
xlim([180 250])
ylim([300 370])
c1 = colorbar;
colormap(cmap)
caxis([0 25])
c1.Label.String = 'Depth (m)';         %Colorbar label

%% - Find Cross-Section Angle

load('\savloc\NOAA_HighResBathy.mat');

%UD Cross
x1 = 216.753;   y1 = 319.825;   x2 = 208.949;   y2 = 331.991;
alongslope = (y2-y1)/(x2-x1);
crossslope = -(1/alongslope);

brandyx = 216.011;
brandyy = 327.993;
uppointx = 221.5;
uppointy = (crossslope*(uppointx-brandyx))+brandyy;
dnpointx = 208.5;
dnpointy = (crossslope*(dnpointx-brandyx))+brandyy;

figure
pcolor(xr,yr,bathy)
%pcolor(xbathy,ybathy,HRbathy)
shading('flat')
hold on
plot(xr(ix,iy),yr(ix,iy),'rp')
plot([x1 x2],[y1 y2],'r-')
plot([dnpointx uppointx],[dnpointy uppointy],'r-')
xlim([180 250])
ylim([300 370])
c1 = colorbar;
%colormap(cmap)
caxis([0 25])
%caxis([-30 0])
c1.Label.String = 'Depth (m)';

%% - Create Line and Find Closest Points

diagpoints = 14;

outxline=linspace(dnpointx,uppointx,diagpoints);
outyline=linspace(dnpointy,uppointy,diagpoints);

for ii = 1:length(outxline)
    r=( (xr-outxline(ii)).^2+(yr-outyline(ii)).^2);
    [ry, iy]= min(r,[],2);
    [rx, ix]= min(ry);
    iy=iy(ix);
    UDxrloc(ii) = xr(ix,iy);    UDyrloc(ii) = yr(ix,iy);
    UDxind(ii) = ix;    UDyind(ii) = iy;
end

figure
pcolor(xr,yr,bathy)
shading('flat')
hold on
plot([dnpointx uppointx],[dnpointy uppointy],'r-')
plot(UDxrloc,UDyrloc,'k-')
xlim([180 250])
ylim([300 370])
c1 = colorbar;
%colormap(cmap)
caxis([0 25])
c1.Label.String = 'Depth (m)';

%% - Download Data from UD Model and Rotate It

if spring == 1
    UDpath = '\savloc\ocean_his_0010.nc';
    ot = nc_varget(UDpath,'ocean_time');
    t0=datenum(2016,9,1,0,0,0);                 %Initial time for the ROMS run (for 2008 - 2008,1,1,0,0,0)
    otc = (ot/3600/24)+t0;
    zeta = nc_varget(UDpath,'zeta');
    timewindtidavg = 509:558;
else
    UDpath = '\savloc\ocean_his_0010.nc';
    ot = nc_varget(UDpath,'ocean_time');
    t0=datenum(2016,9,1,0,0,0);                 %Initial time for the ROMS run (for 2008 - 2008,1,1,0,0,0)
    otc1 = (ot/3600/24)+t0;
    z1 = nc_varget(UDpath,'zeta');
    UDpath = '\savloc\ocean_his_0011.nc';
    ot = nc_varget(UDpath,'ocean_time');
    t0=datenum(2016,9,1,0,0,0);                 %Initial time for the ROMS run (for 2008 - 2008,1,1,0,0,0)
    otc2 = (ot/3600/24)+t0;
    z2 = nc_varget(UDpath,'zeta');
    otc = cat(1,otc1,otc2);
    zeta = cat(1,z1,z2);
    timewindtidavg = 3:52;
end

uUD = nc_varget(UDpath,'u');
uUD = squeeze(uUD(timewindtidavg,:,:,:));
vUD = nc_varget(UDpath,'v');
vUD = squeeze(vUD(timewindtidavg,:,:,:));
gangUD = nc_varget(UDpath,'angle');
if salt == 1
    sUD = nc_varget(UDpath,'salt');
    sUD = squeeze(sUD(timewindtidavg,:,:,:));
end

for ita = 1:size(uUD,1)
    for id = 1:10
        for ia = 1:300
            for ja = 1:148
                urho(ia,ja) = 0.5*(uUD(ita,id,ia,ja)+uUD(ita,id,ia,ja+1));
            end
        end
        fixu = zeros(300,1);
        urhofull = cat(2,fixu,urho,fixu);
        urhofull2U(ita,id,:,:) = urhofull;
    
        for ib = 1:298
            for jb = 1:150
                vrho(ib,jb) = 0.5*(vUD(ita,id,ib,jb)+vUD(ita,id,ib+1,jb));
            end
        end
        fixv = zeros(1,150);
        vrhofull = cat(1,fixv,vrho,fixv);
        vrhofull2U(ita,id,:,:) = vrhofull;
    
        for i = 1:300
            for j = 1:150
                uago(i,j)=(urhofull(i,j)*cos(gangUD(i,j)))-(vrhofull(i,j)*sin(gangUD(i,j)));
                vago(i,j)=(vrhofull(i,j)*cos(gangUD(i,j)))+(urhofull(i,j)*sin(gangUD(i,j)));
            end
        end
    
        uago(1,:)=0; uago(300,:)=0; uago(:,1)=0; uago(:,150)=0;
        vago(1,:)=0; vago(300,:)=0; vago(:,1)=0; vago(:,150)=0;

        urrUD(ita,id,:,:) = uago;
        vrrUD(ita,id,:,:) = vago;
    
    end    
end

uUDg = urhofull2U(:,10,:,:);      vUDg = vrhofull2U(:,10,:,:);
uUDr = squeeze(urrUD(:,10,:,:));  vUDr = squeeze(vrrUD(:,10,:,:));
uUDrf = urrUD;                    vUDrf = vrrUD;
clear urho urhofull urhofull2U vrho vrhofull vrhofull2U uago vago urrUD vrrUD fixu fixv

uUDrt = zeros(length(timewindtidavg),10,length(UDxind));
vUDrt = zeros(length(timewindtidavg),10,length(UDxind));
sUDrt = zeros(length(timewindtidavg),10,length(UDxind));

for it = 1:length(timewindtidavg)
    for id = 1:10
        for i = 1:length(UDxind)
            uUDrt(it,id,i) = uUDrf(it,id,UDxind(i),UDyind(i));
            vUDrt(it,id,i) = vUDrf(it,id,UDxind(i),UDyind(i));
            if salt == 1
                sUDrt(it,id,i) = sUD(it,id,UDxind(i),UDyind(i));
            end
        end
    end
end

%% - Rotate the Cross-Section

invgang = atan2((uppointy-dnpointy),(uppointx-dnpointx));

uacross = zeros(size(uUDrt,1),size(uUDrt,2),size(uUDrt,3));
valong = zeros(size(uUDrt,1),size(uUDrt,2),size(uUDrt,3));

for it2 = 1:size(uUDrt,1)
    for id2 = 1:size(uUDrt,2)
        for iii = 1:size(uUDrt,3)
            uacross(it2,id2,iii)=(uUDrt(it2,id2,iii)*cos(invgang))+(vUDrt(it2,id2,iii)*sin(invgang));
            valong(it2,id2,iii)=(vUDrt(it2,id2,iii)*cos(invgang))-(uUDrt(it2,id2,iii)*sin(invgang));
        end
    end
end

%% - Calculate Distance and Find Depth

bathycross = zeros(1,length(UDxind));
for i2 = 1:length(UDxind)
    bathycross(i2) = bathy(UDxind(i2),UDyind(i2));
end

cd1 = zeros(1,length(UDxind));
for i3 = 1:(length(UDxind)-1)
    x1 = UDxrloc(i3);   x2 = UDxrloc(i3+1);
    y1 = UDyrloc(i3);   y2 = UDyrloc(i3+1);
    cd1(i3) = sqrt(((x2-x1)^2)+((y2-y1)^2));
end
crossdist = cumsum(cd1);

%% - Plot Tidal Average Surface Data

acrosstidavg = mean(uacross,1);     acrosstidavg = squeeze(acrosstidavg(1,:,:));
alongtidavg = mean(valong,1);       alongtidavg = squeeze(alongtidavg(1,:,:));
if salt == 1
    salttidavg = mean(sUDrt,1);        salttidavg = squeeze(salttidavg(1,:,:));
end

figure
set(gcf,'Position',[150   150   1400   600])
subplot(2,1,1)
plot(crossdist,squeeze(alongtidavg(10,:)),'k-')
hold on
plot(crossdist,squeeze(acrosstidavg(10,:)),'r-')
grid on
xlim([0 ceil(max(crossdist))])
ylim([-0.5 0.5])
yticks([-0.5:0.1:0.5])
legend('Along-Channel','Across-Channel')
if spring == 1
    title('UD Model Projected Velocities - Tidal Average - Spring')
else
    title('UD Model Projected Velocities - Tidal Average - Neap')
end
xlabel('Cross-Channel Distance (km)')
ylabel('Velocity (m s^{-1})')

subplot(2,1,2)
plot(crossdist,-bathycross,'k-')
xlim([0 ceil(max(crossdist))])
ylim([-25 0])
grid on
title('Bathymetry')
xlabel('Cross-Channel Distance (km)')
ylabel('Depth (m)')

if isv == 1
    print([dirname,'SurfaceVel'],'-dpng');
end

%% - Download W velocities and Set-Up Vertical Depths

wa = nc_varget(UDpath,'w');
wa = squeeze(wa(timewindtidavg,:,:,:));
for ivi = 1:length(timewindtidavg)
   for id3 = 1:11
        for ivii = 1:length(UDxind)
            wacut(ivi,id3,ivii) = wa(ivi,id3,UDxind(ivii),UDyind(ivii));
        end
    end
end
waavg = mean(wacut,1);            waavg = squeeze(waavg(1,:,:));
for iw = 1:size(waavg,2)
    for i = 1:10
        wa1 = waavg(i,iw);     wa2 = waavg(i+1,iw);
        warray(i,iw)=(wa1+wa2)/2;
    end
end

%SCOORD
 % - Alan Input
    theta_s = 3;   theta_b = 0.4;   Tcline = 0.0;   N=10;  %Values from ROMS
    kgrid = 0;      %Rho grid = 0; w grid = 1;
    h = bathycross;
    
    % - Scoord parameters
    c1=1.0;     c2=2.0;     p5=0.5;     Np=N+1;     ds=1.0/N;
    hmin=min(min(bathy));   hmax=max(max(bathy));   hc=min(hmin,Tcline);
    [Mp,Lp]=size(h);
    
    % - Define S-Curves at vertical RHO- or W-points (-1 < sc < 0).
    if (kgrid == 1)
      Nlev=Np;
      lev=0:N;
      sc=-c1+lev.*ds;
    else
      Nlev=N;
      lev=1:N;
      sc=-c1+(lev-p5).*ds;
    end
    Ptheta=sinh(theta_s.*sc)./sinh(theta_s);
    Rtheta=tanh(theta_s.*(sc+p5))./(c2*tanh(p5*theta_s))-p5;
    Cs=(c1-theta_b).*Ptheta+theta_b.*Rtheta;
    Cd_r=(c1-theta_b).*cosh(theta_s.*sc)./sinh(theta_s)+ ...
         theta_b./tanh(p5*theta_s)./(c2.*(cosh(theta_s.*(sc+p5))).^2);
    Cd_r=Cd_r.*theta_s;
    
    % - Compute depths at requested grid section
    %   Removed zero free-surface assumption!!!
      z=zeros(Lp,Nlev);
      zeta0 = zeros(Lp,Nlev);
      for k=1:Nlev
        z(:,k)=zeta0(:,1).*(c1+sc(k))+hc.*sc(k)+ ...
                 (h(1,:)'-hc).*Cs(k);
      end
      bathyarray = -z';

      crossdistf = ones(size(bathyarray,1),size(bathyarray,2)).*crossdist;

%% - Plot Cross-Section and Check Angle

if anglecheck == 1

    angvary = -30:0.5:30;
    deg2rad = angvary.*(pi/180);
    for i = 1:length(deg2rad)
        invgang2(i,:) = invgang-(deg2rad(i));
    end

    UDcrossu = mean(uUDrt,1);   UDcrossu = squeeze(UDcrossu(1,:,:));
    UDcrossv = mean(vUDrt,1);   UDcrossv = squeeze(UDcrossv(1,:,:));

    uacrossang = zeros(length(angvary),size(UDcrossu,1),size(UDcrossu,2));
    valongang = zeros(length(angvary),size(UDcrossu,1),size(UDcrossu,2));

    for iang = 1:length(angvary)
        for id3 = 1:size(UDcrossu,1)
            for iii = 1:size(UDcrossu,2)
                uacrossang(iang,id3,iii)=(UDcrossu(id3,iii).*cos(invgang2(iang)))+(UDcrossv(id3,iii).*sin(invgang2(iang)));
                valongang(iang,id3,iii)=(UDcrossv(id3,iii).*cos(invgang2(iang)))-(UDcrossu(id3,iii).*sin(invgang2(iang)));
            end
        end
    end

    figure
    set(gcf,'Position',[150   150   1400   400])
    cnt = 1001;
    for ian2 = 1:size(valongang,1)
        pcolor(crossdistf,-bathyarray,squeeze(valongang(ian2,:,:)))
        shading('flat')
        hold on
        conval = -0.5:0.05:0.5;
        contour(crossdistf,-bathyarray,squeeze(valongang(ian2,:,:)),conval,'k-','ShowText','on')
        plot(crossdist,-bathycross,'k-')
        grid on
        if spring == 1
            title(['UD Model Velocities Projected - Tidal Avg - Spring- Ang: ',num2str(angvary(ian2))])
        else
            title(['UD Model Velocities Projected - Tidal Avg - Neap - Ang: ',num2str(angvary(ian2))])
        end
        xlabel('Cross-Channel Distance (km)','FontSize',20)
        ylabel('Depth (m)','FontSize',20)
        ylim([min(-bathycross) 0])
        ax = gca; 
        ax.FontSize = 20;
        c1 = colorbar;
        cmap = cmocean('balance',40);
        colormap(cmap);
        caxis([-0.5 0.5])   %ETM 7
        c1.Label.String = 'Along-Channel Velocity (m s^{-1})';
        asfV = 3;           asfw = 500;
        %quiver(crossdistf,-bathyarray,(DBOFSurecovf.*asfV),(warray.*asfw),'r')
        distarrayqu = cat(1,crossdistf,crossdistf(1,:),crossdistf(1,:));
        bathyarrayqu = cat(1,-bathyarray,(-bathyarray(1,:)-7),(-bathyarray(1,:)-9));
        blankadd = zeros(2,size(uacrossang,3)).*nan;
        uarrayqu = cat(1,squeeze(uacrossang(ian2,:,:)),blankadd);
        warrayqu = cat(1,warray,blankadd);
        uarrayqu(11,13) = 0;        warrayqu(11,13) = 0.0005;
        uarrayqu(12,13) = 0.05;     warrayqu(12,13) = 0.0000;
        quiver(distarrayqu,bathyarrayqu,(uarrayqu.*asfV),(warrayqu.*asfw),'r')
        if isv == 1
            print([dirname,'\AngleCheck\',num2str(cnt)],'-dpng');
            clf
            cnt = cnt+1;
        end
    end
else
    figure
    set(gcf,'Position',[150   150   1400   400])
    if salt == 1
        pcolor(crossdistf,-bathyarray,salttidavg)
    else
        pcolor(crossdistf,-bathyarray,alongtidavg)
    end
    shading('flat')
    hold on
    if salt == 1
        conval = 0:1:36;
        contour(crossdistf,-bathyarray,salttidavg,conval,'k-','ShowText','on')
    else
        conval = -0.5:0.05:0.5;
        contour(crossdistf,-bathyarray,alongtidavg,conval,'k-','ShowText','on')
    end
    plot(crossdist,-bathycross,'k-')
    grid on
%     if spring == 1
%         title('UD Model Project Velocities - Tidal Avg - Spring')
%     else
%         title('UD Model Project Velocities - Tidal Avg - Neap')
%     end
    xlabel('Cross-Channel Distance (km)','FontSize',20)
    ylabel('Depth (m)','FontSize',20)
    ylim([min(-bathycross) 0])
    ax = gca; 
    ax.FontSize = 20;
    c1 = colorbar;
    if salt == 1
        cmap = cmocean('haline');
        %cmap = cmocean('balance',40);
        colormap(cmap);
        caxis([22 32])   %ETM 7
        c1.Label.String = 'Salinity (PSU)';
    else
        cmap = cmocean('balance',40);
        colormap(cmap);
        caxis([-0.5 0.5])   %ETM 7
        c1.Label.String = 'Along-Channel Velocity (m s^{-1})';
    end
    asfV = 1;           asfw = 500;
    %quiver(crossdistf,-bathyarray,(DBOFSurecovf.*asfV),(warray.*asfw),'r')
    distarrayqu = cat(1,crossdistf,crossdistf(1,:),crossdistf(1,:));
    bathyarrayqu = cat(1,-bathyarray,(-bathyarray(1,:)-7),(-bathyarray(1,:)-9));
    blankadd = zeros(2,size(acrosstidavg,2)).*nan;
    uarrayqu = cat(1,acrosstidavg,blankadd);
    warrayqu = cat(1,warray,blankadd);
    uarrayqu(11,9) = 0;        warrayqu(11,9) = 0.0005;
    uarrayqu(12,9) = 0.05;     warrayqu(12,9) = 0.0000;
    %quiver(distarrayqu,bathyarrayqu,(uarrayqu.*asfV),(warrayqu.*asfw),'r','LineWidth',2)
    quiver(distarrayqu,bathyarrayqu,(uarrayqu.*asfV),(warrayqu.*asfw),'w','LineWidth',4)
    quiver(distarrayqu,bathyarrayqu,(uarrayqu.*asfV),(warrayqu.*asfw),'r','LineWidth',1.75)
    
    if isv == 1
        if salt == 1
            print([dirname,'FullCrossSection_Salinity'],'-dpng');
        else
            print([dirname,'FullCrossSection_AlongVel'],'-dpng');
        end
    end
end

if isv == 1
    %Save Tidally Varying Data
    snvary = [dirname,'TidallyVaryingData'];    
    save(snvary,'crossdist','bathyarray','uacross','valong')

end

toc