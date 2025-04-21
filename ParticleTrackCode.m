%%%Creates particles and tracks their movement in the Delaware Bay ROMS flow field
%%%Update - Corrects for curvilinear
%------------------------------------------------------------------------
%Author: R Alan Mason --> Based on original code by T. Kukulka
%Created: 7/6/21 
%Edited: 12/21/21
%------------------------------------------------------------------------
%***Last edit: Particle field Initialization

%%%Flags and Parameters----------------------------------------------------
    %Sets a flag for what month is being used------------------------------
    IdealRST=1;             % If using the rst file that begins with Mar 1
    % -- misc parameters --------------------------------------------------
    dt = 1*300;             % particle tracking time step in s(!) ---> 300 seconds = 5 min
    numdays = 10;
    di=15;                   % used for minimum search to determine (ix,iy) ---> Original = 5
    bounce=0;               % apply bounce if =1
    isv=1;                  % save image files if isv==1 for movie, use with mk_movie.m
    plotsamplet = 0;        % Have sample stations highlight when sampled?
    dirname = '\mv\';        % save image files here
    cnt=1000;               % count of save image files
    %stepsi = 0.1;          % In kilometers - Vaild Values = 0.025, 0.05, 0.1, 0.25, 0.5, 1.0
    stepsi = 0.25;
    hashtagbinning=1;       % Flag that turns binning off (0) and on (1)
    binsize = stepsi*3;     % Sets the bin size (in square km)
    binshow = 0;            % Flag that turns showing the bins off (0) and on (1) ---> Off is recommended for small bin sizes!!!
    sto = 0;                % Flag that turns Stochastic Forcing off (0) and on (1)
    edD = (4.2/(1000^2));   % Sets Eddy Diffusivity (km scale) - This set (m2/s over km2 to correct units) %%Wind=0.0148 m2/s
    smallsave = 1;          % Flag that sets up saving the xs/ys files in chunks
    lensave = 12*(3600/dt); % Time in hours per file converted to number of time steps
    smallname = '\PartTrack\'; % Name for saved file
    usenormbin = 0;         % Flag that allows for normalization of histogram bin values
    binvalmax = 3660;      % Maximum binvalue for 1.3 million particles 
    mymap=[0.28627 0.28627 0.30196; 0.52941 0.55294 0.57255; 0.89020 0.89020 0.80392; 0.87451 0.40392 0.27843; 0.70588 0.18431 0.19608];
    suvname='C:\Users\ramason\Documents\MATLAB\DE_Bay\suvfiles\IdealizedTests\Riv_Avg_Wind_0\suv';
    %Particles Parameters--------------------------------------------------
    addpart = 0;            % Flag that adds particles to the code
    addstep = 3;            % Hours between adding particles
    addelay = 0;            % Hours to delay before adding particles
    
    upbayxline = linspace(170,177,22);
    upbayyline = linspace(396,400,22);
    [xnewup,ynewup] = meshgrid(upbayxline,upbayyline);
    [ny,nx]=size(xnewup);
    xnewup=reshape(xnewup,[1,nx*ny]);
    ynewup=reshape(ynewup,[1,nx*ny]);
    xq = [172.8,176.4,174.5,170.9];
    yq  = [396.4,398.9,399.8,397.3];
    [in,on] = inpolygon(xnewup,ynewup,xq,yq);
    xnewup = xnewup(in);
    ynewup = ynewup(in);

    xnew = xnewup;
    ynew = ynewup;
    
    numadd = length(xnew);  % The number of random particles added per timestep
    includeamprand = 0;     % Flag that adds some randomness to the added particle
    valamprand = 0.4;       % Random Amplitude value
    plastic = 0;
    cmap=colormap;
    upside_down_cmap=flipud(cmap);
    cmapX = cool(60); cmapY = spring(30); cmapZ = flipud(parula(15));
    cmapsalt = [cmapX;cmapY;cmapZ];
    %Partial Particle Field------------------------------------------------
    partialpart = 0;        % Flag that allows for the partitioning of the particle field
    partialpart1bound = 1;        % Flag that allows for the partitioning of the particle field
    partshelf = 0;          % Flag that initializes particles on the shelf off the Delaware Coast
    partRho = 0;            % Flag that initializes particles are rho points
    partDrift = 0;
    ToddNode = 0;
    
    if partialpart == 1
        xbUPmin = 184.5;        % Sets lower x-bound of particle field %Central bay = 205, DE = 212, NJ = 232, Mouth = 225, Upper = 184.5 
        xbUPmax = 189.5;        % Sets upper x-bound of particle field %Central bay = 210, DE = 217, NJ = 237, Mouth = 230, Upper = 189.5 
        ybUPmin = 380;        % Sets lower y-bound of particle field %Central bay = 340, DE = 309, NJ = 353, Mouth = 315, Upper = 380
        ybUPmax = 385;        % Sets upper y-bound of particle field %Central bay = 345, DE = 314, NJ = 358, Mouth = 320, Upper = 385

        xbCBmin = 216; %205       % Sets lower x-bound of particle field for the Central bay
        xbCBmax = 231; %210       % Sets upper x-bound of particle field for the Central bay 
        ybCBmin = 333; %340       % Sets lower y-bound of particle field for the Central bay
        ybCBmax = 360; %345       % Sets upper y-bound of particle field for the Central bay

        xbDEmin = 202; %214       % Sets lower x-bound of particle field near Cape Henlopen 
        xbDEmax = 215; %219       % Sets upper x-bound of particle field near Cape Henlopen 
        ybDEmin = 322; %311       % Sets lower y-bound of particle field near Cape Henlopen
        ybDEmax = 342; %316       % Sets upper y-bound of particle field near Cape Henlopen

        xbNJmin = 232;        % Sets lower x-bound of particle field for the NJ shore 
        xbNJmax = 237;        % Sets upper x-bound of particle field for the NJ shore
        ybNJmin = 340;        % Sets lower y-bound of particle field for the NJ shore
        ybNJmax = 358;        % Sets upper y-bound of particle field for the NJ shore

        xbMOmin = 210;        % Sets lower x-bound of particle field for the central Bay Mouth 
        xbMOmax = 235;        % Sets upper x-bound of particle field for the central Bay Mouth 
        ybMOmin = 311;        % Sets lower y-bound of particle field for the central Bay Mouth
        ybMOmax = 321;        % Sets upper y-bound of particle field for the central Bay Mouth
        
        xbDOVmin = 191;        % Sets lower x-bound of particle field near Dover, DE 
        xbDOVmax = 209;        % Sets upper x-bound of particle field near Dover, DE  
        ybDOVmin = 343;        % Sets lower y-bound of particle field near Dover, DE 
        ybDOVmax = 365;        % Sets upper y-bound of particle field near Dover, DE 
        
    elseif partialpart1bound == 1
        %xboundmin = 184; %188;        % Sets lower x-bound of particle field %Central bay = 205, DE = 212, NJ = 232, Mouth = 225, Upper = 184.5 
        %xboundmax = 186; %231;        % Sets upper x-bound of particle field %Central bay = 210, DE = 217, NJ = 237, Mouth = 230, Upper = 189.5 
        %yboundmin = 382; %310;        % Sets lower y-bound of particle field %Central bay = 340, DE = 309, NJ = 353, Mouth = 315, Upper = 380
        %yboundmax = 384; %377;        % Sets upper y-bound of particle field %Central bay = 345, DE = 314, NJ = 358, Mouth = 320, Upper = 385
        %DBOFS Bounds
        xboundmin = 187; %188;        % Sets lower x-bound of particle field %Central bay = 205, DE = 212, NJ = 232, Mouth = 225, Upper = 184.5 
        xboundmax = 236.1; %231;        % Sets upper x-bound of particle field %Central bay = 210, DE = 217, NJ = 237, Mouth = 230, Upper = 189.5 
        yboundmin = 305.3; %310;        % Sets lower y-bound of particle field %Central bay = 340, DE = 309, NJ = 353, Mouth = 315, Upper = 380
        yboundmax = 380;
        
    elseif partshelf == 1
        xboundmin = 235;        % Sets lower x-bound of particle field %Central bay = 205, DE = 212, NJ = 232, Mouth = 225, Upper = 184.5 
        xboundmax = 240;        % Sets upper x-bound of particle field %Central bay = 210, DE = 217, NJ = 237, Mouth = 230, Upper = 189.5 
        yboundmin = 275;        % Sets lower y-bound of particle field %Central bay = 340, DE = 309, NJ = 353, Mouth = 315, Upper = 380
        yboundmax = 280;        % Sets upper y-bound of particle field
    
    end
    
%--------------------------------------------------------------------------

%%%Load in Constant Variables----------------------------------------------
    % Note - These files are a combination of xy_ruv_bath.mat, Coasts.mat, and the excel particle location files
if stepsi == 0.025
    load('C:\Users\ramason\Documents\GitHub\ParticleTrackFiles\ParticleTrack_13_7v2_Initialfiles\StartUpFile_Step25m.mat')
elseif stepsi == 0.05
    load('C:\Users\ramason\Documents\GitHub\ParticleTrackFiles\ParticleTrack_13_7v2_Initialfiles\StartUpFile_Step50m.mat')
elseif stepsi == 0.1
    load('C:\Users\ramason\Documents\GitHub\ParticleTrackFiles\ParticleTrack_13_7v2_Initialfiles\StartUpFile_Step100m.mat')
elseif stepsi == 0.25
    load('C:\Users\ramason\Documents\GitHub\ParticleTrackFiles\ParticleTrack_13_7v2_Initialfiles\StartUpFile_Step250m.mat')
elseif stepsi == 0.5
    load('C:\Users\ramason\Documents\GitHub\ParticleTrackFiles\ParticleTrack_13_7v2_Initialfiles\StartUpFile_Step500m.mat')
else
    load('C:\Users\ramason\Documents\GitHub\ParticleTrackFiles\ParticleTrack_13_7v2_Initialfiles\StartUpFile_Step1000m.mat')
end
%--------------------------------------------------------------------------

%%%Creation of time for ROMS and Particle Tracking-------------------------
    %ROMS time and file output parameters ------------------------------------
    if IdealRST == 1
        t0_month=datenum(2017,3,31,11,0,0);  % ROMS time is "seconds since 2008-01-01 00:00:00" - May 2008
        %dateofmonth = 7 - 1;               % Sets starting date within the month (date-1 to account for the first day being day zero)
        %hourofday = (12/24);                 % Sets starting hour on the selected day above (0-23) %%6
        %dateofmonth = 17 - 1;      %Neap Mid-Month
        %hourofday = (7/24);
        dateofmonth = 23 - 1;       %Spring End of Month - 10 Day
        hourofday = (23/24);
        %dateofmonth = 26 - 1;       %Spring End of Month
        %hourofday = (2/24);
        t0_ROMS=datenum(2017,3,1,0,0,0);
        time = nc_varget('G:\ram_Research\Model_Runs\ramason_runs\Idealized_Tests\Riv_Avg_Wind_0\BayMouthCross\ocean_his_0002.nc','ocean_time');
    else
        t0_month=datenum(2016,9,1,0,0,0);   % ROMS time is "seconds since 2008-01-01 00:00:00"
        dateofmonth = 1 - 1;                % Sets starting date within the month (date-1 to account for the first day being day zero)
        hourofday = (0/24);                 % Sets starting hour on the selected day above (0-23)
        t0_ROMS=datenum(2016,9,1,0,0,0);   % ROMS time is "seconds since 2008-01-01 00:00:00"
    end
    dtime = 3600;                           % Number of seconds for ROMS timestep
    %tlh = 730*ones(13,1); tlh(13) = 24;     % Length of ROMS files in hours
    tlh = 730*ones(25,1); tlh(25) = 24;     % Length of ROMS files in hours
    tlh = cumsum(tlh);                      % List of last times in hours from t0 for each file i=1...13
    %tlh(2:end)=tlh(2:end)-1;
    itimeold=-1;                            % Sets up parameter for later use in the code (tells code to advance hour in ROMS)
    % particle tracking time parameters ---------------------------------------
    nt= numdays*24*round(3600/dt);          % number of particle tracking time steps
    t=dt*(0:nt-1);                          % particle tracking time: Note this is different from ROMS!
    t=t/24/3600+t0_month;                   % in DAYS from the start of the month!
    t=t+dateofmonth+hourofday;              % Adjusts t to start at a specified time in the middle fo the month
    startdatetime = datestr((t(1)));        % Lists the initial time for use in small save
    % set-up addparticle time array
    addadjust = (addelay*round(3600/dt))+1;
    addparthour = addadjust:(addstep*round(3600/dt)):nt;
    %addadjust = (addelay/addstep)+1;
    %addparthour = addparthour(addadjust:end);
    
%--------------------------------------------------------------------------

%%%Load the first ROMS File of the study period----------------------------
th=floor((t(1)-t0_ROMS)*24);                % ROMS time in hours from t0
ithf=find(tlh>=th,1,'first');               % Finds first intance that tlh >= th to place the month the hour occurs in
%ithf=4;
load([suvname,num2str(ithf)]) % Loads in 's','u','v','time','readme','fname' for the initial month
ubart=u; vbart=v; sbart=s; clear u v s      % Renames variables u, v, and s then deletes them
gang = nc_varget('\dataloc\ocean_his_0008.nc','angle');

for ita = 1:size(ubart,1)
    for ia = 1:300
        for ja = 1:148
            urho(ia,ja) = 0.5*(ubart(ita,ia,ja)+ubart(ita,ia,ja+1));
        end
    end
    fixu = zeros(300,1);
    urhofull = cat(2,fixu,urho,fixu);
    for ib = 1:298
        for jb = 1:150
            vrho(ib,jb) = 0.5*(vbart(ita,ib,jb)+vbart(ita,ib+1,jb));
        end
    end
    fixv = zeros(1,150);
    vrhofull = cat(1,fixv,vrho,fixv);

    for i = 1:300
        for j = 1:150
            uago(i,j)=(urhofull(i,j)*cos(gang(i,j)))-(vrhofull(i,j)*sin(gang(i,j)));
            vago(i,j)=(vrhofull(i,j)*cos(gang(i,j)))+(urhofull(i,j)*sin(gang(i,j)));
        end
    end

    uago(1,:)=0; uago(300,:)=0; uago(:,1)=0; uago(:,150)=0;
    vago(1,:)=0; vago(300,:)=0; vago(:,1)=0; vago(:,150)=0;

    urr(ita,:,:) = uago;
    vrr(ita,:,:) = vago;
    
end

ubart = urr; vbart = vrr;
clear urho urhofull vrho vrhofull uago vago urr vrr fixu fixv
%--------------------------------------------------------------------------

%%%Initialize particle locations and integration---------------------------

    %Uses the loaded excel values to generate all particle locations-------
    
    bathy(250:300,1:60)=nan;
    bathynew = zeros(300,150).*NaN;
    for i1 = 3:299
        for i2 = 1:150
            bathymid = bathy(i1,i2);
            bathyup = bathy((i1+1),i2);
            bathydwn = bathy((i1-2),i2);

            if isnan(bathymid) || isnan(bathyup) || isnan(bathydwn)
                bathynew(i1,i2) = NaN;
            else
                bathynew(i1,i2) = bathymid;
            end

        end
    end
    for i1 = 1:300
        for i2 = 2:149
            bathymid = bathy(i1,i2);
            bathyup = bathy(i1,(i2+1));
            bathydwn = bathy(i1,(i2-1));

            if isnan(bathymid) || isnan(bathyup) || isnan(bathydwn)
                bathynew(i1,i2) = NaN;
            else
                %bathynew(i1,i2) = bathymid;
            end

        end
    end

    xlmin = 180; xlmax = 244; xstep = (xlmax-xlmin)/stepsi;
    ylmin = 302; ylmax = 396; ystep = (ylmax-ylmin)/stepsi;
    fullbayxline = linspace(xlmin,xlmax,xstep);
    fullbayyline = linspace(ylmin,ylmax,ystep);
    [xnewup,ynewup] = meshgrid(fullbayxline,fullbayyline);
    [ny,nx]=size(xnewup);
    xfull=reshape(xnewup,[1,nx*ny]);
    yfull=reshape(ynewup,[1,nx*ny]);

    np = length(xfull);

    ixs = zeros(1,np); iys=zeros(1,np);     % Grid Location Array
    clearar = zeros(1,np);
    for ipi=1:np
        x=xfull(1,ipi); y=yfull(1,ipi);
        r=( (xr-x).^2+(yr-y).^2);
        [ry, iy]= min(r,[],2);
        [rx, ix]= min(ry);
        iy=iy(ix);
        ixs(ipi)=ix; iys(ipi)=iy;
        
        bathyval = bathynew(ix,iy);
        if isnan(bathyval)
            clearar(1,ipi) = 1;
        else
            clearar(1,ipi) = 0;
        end
    end

    xfull(clearar == 1) = [];
    yfull(clearar == 1) = [];
    
    clear clearar

   if partialpart == 1
        partcount = 1;
        %%%Set up area arrays
            xpUP = 0; ypUP = 0;
            xpCB = 0; ypCB = 0;
            xpDE = 0; ypDE = 0;
            xpNJ = 0; ypNJ = 0;
            xpMO = 0; ypMO = 0;
            xpDOV = 0; ypDOV = 0;
            %xpEI = 0; ypEI = 0;
        for partclean = 1:length(xfull)
            xclean = xfull(partclean);
            yclean = yfull(partclean);
            if xclean >= xbUPmin && xclean <= xbUPmax && yclean >= ybUPmin && yclean <= ybUPmax
                %if yclean >= ybUPmin && yclean <= ybUPmax
                    xpartial(partcount) = xfull(partclean);
                    ypartial(partcount) = yfull(partclean);
                    xpUP = cat(2,xpUP,partcount);
                    ypUP = cat(2,ypUP,partcount);
                    partcount = partcount+1;
                %end
            elseif xclean >= xbCBmin && xclean <= xbCBmax && yclean >= ybCBmin && yclean <= ybCBmax
                %if yclean >= ybCBmin && yclean <= ybCBmax
                    xpartial(partcount) = xfull(partclean);
                    ypartial(partcount) = yfull(partclean);
                    xpCB = cat(2,xpCB,partcount);
                    ypCB = cat(2,ypCB,partcount);
                    partcount = partcount+1;
                %end
            elseif xclean >= xbDEmin && xclean <= xbDEmax && yclean >= ybDEmin && yclean <= ybDEmax
                %if yclean >= ybDEmin && yclean <= ybDEmax
                    xpartial(partcount) = xfull(partclean);
                    ypartial(partcount) = yfull(partclean);
                    xpDE = cat(2,xpDE,partcount);
                    ypDE = cat(2,ypDE,partcount);
                    partcount = partcount+1;
                %end
            elseif xclean >= xbNJmin && xclean <= xbNJmax && yclean >= ybNJmin && yclean <= ybNJmax
                %if yclean >= ybNJmin && yclean <= ybNJmax
                    xpartial(partcount) = xfull(partclean);
                    ypartial(partcount) = yfull(partclean);
                    xpNJ = cat(2,xpNJ,partcount);
                    ypNJ = cat(2,ypNJ,partcount);
                    partcount = partcount+1;
                %end
            elseif xclean >= xbMOmin && xclean <= xbMOmax && yclean >= ybMOmin && yclean <= ybMOmax
                %if yclean >= ybMOmin && yclean <= ybMOmax
                    xpartial(partcount) = xfull(partclean);
                    ypartial(partcount) = yfull(partclean);
                    xpMO = cat(2,xpMO,partcount);
                    ypMO = cat(2,ypMO,partcount);
                    partcount = partcount+1;
                %end
            elseif xclean >= xbDOVmin && xclean <= xbDOVmax && yclean >= ybDOVmin && yclean <= ybDOVmax
                %if yclean >= ybDOVmin && yclean <= ybDOVmax
                    xpartial(partcount) = xfull(partclean);
                    ypartial(partcount) = yfull(partclean);
                    xpDOV = cat(2,xpDOV,partcount);
                    ypDOV = cat(2,ypDOV,partcount);
                    partcount = partcount+1;
                %end
            %elseif xclean >= xbEImin && xclean <= xbEImax && yclean >= ybEImin && yclean <= ybEImax
                %if yclean >= ybEImin && yclean <= ybEImax
                    %xpartial(partcount) = xfull(partclean);
                    %ypartial(partcount) = yfull(partclean);
                    %xpEI = cat(2,xpEI,partcount);
                    %ypEI = cat(2,ypEI,partcount);
                    %partcount = partcount+1;
                %end
            end
        end
        xfull = xpartial;
        yfull = ypartial;
        %%%Clean up area arrays
            xpUP = xpUP(2:end); ypUP = ypUP(2:end);
            xpCB = xpCB(2:end); ypCB = ypCB(2:end);
            xpDE = xpDE(2:end); ypDE = ypDE(2:end);
            xpNJ = xpNJ(2:end); ypNJ = ypNJ(2:end);
            xpMO = xpMO(2:end); ypMO = ypMO(2:end);
            xpDOV = xpDOV(2:end); ypDOV = ypDOV(2:end);
            %xpEI = xpEI(2:end); ypEI = ypEI(2:end);
    elseif partialpart1bound == 1
        for partclean = 1:length(xfull)
            xclean = xfull(partclean);
            yclean = yfull(partclean);
            if xclean < xboundmin || xclean > xboundmax
                xfull(partclean) = 0;
                yfull(partclean) = 0;
            elseif yclean < yboundmin || yclean > yboundmax
                xfull(partclean) = 0;
                yfull(partclean) = 0;  
            end
        end
        xfull(xfull==0)=[];
        yfull(yfull==0)=[];
        
        %xfull = xfull+5;
        %yfull = yfull-75;
    else
         %Nothing to be done
    end
    
    %Creates parameters and matrices that are used for full tracks---------
    np = length(yfull);         % Number of Particles
    npold = np;                 % Set initial value - Used when addpart == 1
    nporiginal = np;            % Saves initial number of particles
    if smallsave == 1
        ys = zeros(lensave,np); % Y-Array to house all particle positions through time
        ys(1,:) = yfull;        % Makes row 1 the initial positions
        xs = zeros(lensave,np); % Same as ys
        xs(1,:) = xfull;
    else
        ys = zeros(nt,np);      % Y-Array to house all particle positions through time
        ys(1,:) = yfull;        % Makes row 1 the initial positions
        xs = zeros(nt,np);      % Same as ys
        xs(1,:) = xfull;
    end
    upre = zeros(np,1);         % Sets up parameters for later use in the code (Previous u/v value used in Adams-Bashforth)
    vpre = zeros(np,1);     
    savnum = 0;                 % Sets up value to be used in numbering saved files
    
    savcorrect = 0;             % Sets up value to be used to correct it value for subsetted xs/ys arrays
    %Finds the closest grid point (ix,iy)----------------------------------
    ixs = zeros(1,np); iys=zeros(1,np);     % Grid Location Array
    for ip=1:np
        x=xs(1,ip); y=ys(1,ip);
        r=( (xr-x).^2+(yr-y).^2);
        [ry, iy]= min(r,[],2);
        [rx, ix]= min(ry);
        iy=iy(ix);
        ixs(ip)=ix; iys(ip)=iy;
    end
%--------------------------------------------------------------------------
    
%%%Sets up plotting parameters---------------------------------------------
fig1=figure;                                    % Used for plotting main figure
set(gcf,'PaperSize',[13.33/2 7.5])              % Sets scale of figure so it does not become distorted
set(gcf,'Color','k')                            % Sets plot background color to be black
set(gcf, 'InvertHardcopy', 'off')               % Tells Matlab to save the figure background color and not to change it to white
if hashtagbinning == 1
    set(gcf,'Position',[378   302   950   653]) % Sets size of figure for plotting both Salinity and Tile Histogram
else
    set(gcf,'Position',[378   302   400   653]) % Sets size of figure for only plotting Salinity
end
fig2=figure;                                    % Used to plot intermediate histogram plots
        % Side Note - For making movies SEE Particle_Track_video.m!!
%--------------------------------------------------------------------------

%%%Code for setting up historgram substitute-------------------------------
xr5=[260,260,260;289,289,289;318,318,318;347,347,347;376,376,376;405,405,405];
[xxrold, yxrold] = meshgrid(1:size(xr5,2),1:size(xr5,1));
[xxrq, yxrq] = meshgrid(linspace(1, size(xr5,2), 566),linspace(1,size(xr5,1),966));
bigxr6=interp2(xxrold, yxrold,xr5,xxrq, yxrq);
x6=bigxr6';

yr5=[180,197,214,231,248,265;180,197,214,231,248,265;180,197,214,231,248,265];
[xyrold, yyrold] = meshgrid(1:size(yr5,2),1:size(yr5,1));
[xyrq, yyrq] = meshgrid(linspace(1, size(yr5,2), 566),linspace(1,size(yr5,1),966));
bigyr6=interp2(xyrold, yyrold,yr5,xyrq, yyrq);
y6=bigyr6';

%--------------------------------------------------------------------------
%%
tic
% ------- time integration loop -------------------------------------------
for it=1:nt-1
    
    %-- get data from ROMS ------------------------------------------------
    ithfold=ithf;
    % first interpolate all points over time, seems redundant but good if many particles (locations) are used and if spatial interpolation expensive
    th = (t(it)-t0_ROMS)*24;            % current time in h from t0
        Ndecimal = 7;
        f = 10.^Ndecimal; 
        th = round(f*th)/f;
    tt=th*3600;                         % current time in s from t0 ---> Used in ROMS time weight
    time1h = floor(th);                 % time1h<=th<=time2h, Time in hours from t0
    time1 = time1h*3600;                % time1<=tt<=time2,   Time in sec from t0, done to be consistent with ROMS time(itime) ---> Used in ROMS time weight
    ithf=find(tlh>=time1h,1,'first');   % file number
    itime=time1h+1; if ithf>1, itime=time1h-tlh(ithf-1); end
    if itimeold==itime
        % keep same ubar1 and ubar2 as for last it, nothing to do
    elseif itime==length(time)
        ithfnew = ithfold+1;
        ubar1=squeeze(ubart(length(time),:,:));
        vbar1=squeeze(vbart(length(time),:,:));
        sbar1=squeeze(sbart(length(time),:,:));
        timeold=time;
        warning('load new file, variable ''time'' not correct for next dtROMS');
        load([suvname,num2str(ithfnew)])
        ubart=u; vbart=v; sbart=s; clear u v s
        for ita = 1:size(ubart,1)
            for ia = 1:300
                for ja = 1:148
                    urho(ia,ja) = 0.5*(ubart(ita,ia,ja)+ubart(ita,ia,ja+1));
                end
            end
            fixu = zeros(300,1);
            urhofull = cat(2,fixu,urho,fixu);
            for ib = 1:298
                for jb = 1:150
                    vrho(ib,jb) = 0.5*(vbart(ita,ib,jb)+vbart(ita,ib+1,jb));
                end
            end
            fixv = zeros(1,150);
            vrhofull = cat(1,fixv,vrho,fixv);

            for i = 1:300
                for j = 1:150
                    uago(i,j)=(urhofull(i,j)*cos(gang(i,j)))-(vrhofull(i,j)*sin(gang(i,j)));
                    vago(i,j)=(vrhofull(i,j)*cos(gang(i,j)))+(urhofull(i,j)*sin(gang(i,j)));
                end
            end

            uago(1,:)=0; uago(300,:)=0; uago(:,1)=0; uago(:,150)=0;
            vago(1,:)=0; vago(300,:)=0; vago(:,1)=0; vago(:,150)=0;

            urr(ita,:,:) = uago;
            vrr(ita,:,:) = vago;

        end

        ubart = urr; vbart = vrr;
        clear urho urhofull vrho vrhofull uago vago urr vrr fixu fixv
        
        ubar2=squeeze(ubart(1,:,:));
        vbar2=squeeze(vbart(1,:,:));
        sbar2=squeeze(sbart(1,:,:));
        [num2str(ithfold),'  ',num2str(ithfnew),'  ',num2str(itime),'  ',datestr(t(it)),'  ',num2str(it)]
    else
        ubar1=squeeze(ubart(itime,:,:));
        vbar1=squeeze(vbart(itime,:,:));
        sbar1=squeeze(sbart(itime,:,:));
        ubar2=squeeze(ubart(itime+1,:,:));
        vbar2=squeeze(vbart(itime+1,:,:));
        sbar2=squeeze(sbart(itime+1,:,:));
    end
    itimeold=itime;
    wt = (tt-time1)/dtime;              % Time weight for ROMS file
    ubar = (1-wt)*ubar1 + wt*ubar2;     % Weighs current and next ubar values based on proxmity to hour mark
    vbar = (1-wt)*vbar1 + wt*vbar2;     % Weighs current and next ubar values based on proxmity to hour mark
    sbar = (1-wt)*sbar1 + wt*sbar2;     % Weighs current and next ubar values based on proxmity to hour mark
    
    %%%Code that adds in new particles while code is running
    if addpart ==1 && ismember(it,addparthour)
        if smallsave == 1
            newcolx = zeros(lensave,numadd); newcoly = zeros(lensave,numadd);     % Creates columns for new particles
        else
            newcolx = zeros(nt,numadd); newcoly = zeros(nt,numadd);
        end
        for numaddpart = 1:numadd
            if includeamprand == 1
                newcolx((it-savcorrect),numaddpart) = xnew(numaddpart)+(valamprand*randn);   % Sets initial position for the particle at the current time --> Includes slight randomness in position
                newcoly((it-savcorrect),numaddpart) = ynew(numaddpart)+(valamprand*randn);
            else
                newcolx((it-savcorrect),numaddpart) = xnew(numaddpart);                      % Sets initial position for the particle at the current time --> No randomness in position
                newcoly((it-savcorrect),numaddpart) = ynew(numaddpart);
            end
            r=( (xr-newcolx((it-savcorrect),numaddpart)).^2+(yr-newcoly((it-savcorrect),numaddpart)).^2); %Finds closes grid position and adds it to the ixs/iys array
            [ry, iy]= min(r,[],2);
            [rx, ix]= min(ry);
            iy=iy(ix);
            ixs(np+numaddpart)=ix; iys(np+numaddpart)=iy;
        end
        xs = [xs newcolx]; ys = [ys newcoly];       % Adds new particle columns to the xs/ys array
        npold = np;                                 % Saves the old np value for later use
        np = length(squeeze(xs(1,:)));              % Updates the np value so all particles are evaluated in the next step
    end
    
    for ip=1:np
        
        %-- compute x-> x+dx ----------------------------------------------
        x=xs((it-savcorrect),ip); y=ys((it-savcorrect),ip);   %current position of particle ip   
        ix=ixs(ip); iy=iys(ip);     %approximate closest grid point indices
        
        %-- get u, v at (x,y) for each particle ---------------------------
        % start with u
        % first get nearest grid point (xu(ix,iy),yu(ix,iy))
        ix1 = ix; iy1=iy;   % pick a square with (2 di + 1)^2 grid points centered at (ix1,iy1)
        iix = ix1-di:ix1+di; iiy = iy1-di:iy1+di;
        ncutx=iix<1; ncuty=iiy<1;
        iix(iix>nxu|ncutx)=[]; iiy(iiy>nyv|ncuty)=[];   % search over those grid points --> Find values inside the domain
        r=( (xr(iix,iiy)-x).^2+(yr(iix,iiy)-y).^2);     % note: r is computed on subdomain
        [ry, iy]= min(r,[],2); 
        [rx, ix]= min(ry); iy=iy(ix);
        ix=ix+ix1-1-di; iy=iy+iy1-1-di;     % convert indices from small square back to whole domain
        ix=ix+sum(ncutx); iy=iy+sum(ncuty); % correct for r cut at small indices
        ixi=[ix,ix+1,ix-1,ix,ix]; iyi=[iy,iy,iy,iy+1,iy-1]; % neighboring grid points
        %ixi=[ix,ix+1,ix-1,ix,ix+1,ix-1,ix,ix+1,ix-1,]; iyi=[iy,iy,iy,iy+1,iy+1,iy+1,iy-1,iy-1,iy-1];
        ii = ixi>nxu | iyi>nyv | ixi<1 | iyi<1; % remove outside grid points
        ixi(ii)=[]; iyi(ii)=[];
        tmp = diag(ubar(ixi,iyi));
        ii = isnan(tmp);  % ii(1) is ubar closest to (x,y)
        if ii(1), bounce=1; bouncu=1; iib=ii; ixb=ixi;  iyb=iyi; end
        ixi(ii)=[]; iyi(ii)=[]; % outside domain- need to make sure particle not on land!
        % weight flow linearly based on distances
        nixi=length(ixi);
        rs=zeros(nixi,1);
        w=zeros(nixi,1);
        for i=1:nixi
            ixr = ixi(i)-ix1+1+di-sum(ncutx);   % convert to local r index
            iyr = iyi(i)-iy1+1+di-sum(ncuty);
            rs(i)=r(ixr,iyr);
            %w(i) = prod(rs( i~=1:nixi) );      %Edit - RAM 1/15/21
        end
        for i=1:nixi
            w(i) = prod(rs( i~=1:nixi) );       %Edit - RAM 1/15/21
        end
        w = w/sum(w);
        u=0;
        for i=1:nixi
            u=u+w(i)*ubar(ixi(i),iyi(i));
        end 
        % do same for v
        ix1 = ix; iy1=iy;   % pick a square with (2 di + 1)^2 grid points centered at (ix1,iy1)
        iix = ix1-di:ix1+di; iiy = iy1-di:iy1+di;
        ncutx=iix<1; ncuty=iiy<1;
        iix(iix>nxu|ncutx)=[]; iiy(iiy>nyv|ncuty)=[];   % search over those grid points
        r=( (xr(iix,iiy)-x).^2+(yr(iix,iiy)-y).^2);     % note: r is computed on subdomain
        [ry, iy]= min(r,[],2); [rx, ix]= min(ry); iy=iy(ix);
        ix=ix+ix1-1-di; iy=iy+iy1-1-di;     % convert indices from small square back to whole domain
        ix=ix+sum(ncutx); iy=iy+sum(ncuty); % correct for r cut at small indices
        ixi=[ix,ix+1,ix-1,ix,ix]; iyi=[iy,iy,iy,iy+1,iy-1];
        %ixi=[ix,ix+1,ix-1,ix,ix+1,ix-1,ix,ix+1,ix-1,]; iyi=[iy,iy,iy,iy+1,iy+1,iy+1,iy-1,iy-1,iy-1];
        ii = ixi>nxu | iyi>nyv | ixi<1 | iyi<1; % remove outsie grid points
        ixi(ii)=[]; iyi(ii)=[];
        tmp = diag(vbar(ixi,iyi));
        ii = isnan(tmp);  % ii(1) is vbar closest to (x,y)
        if ii(1), bounce=1; bouncu=0; iib=ii; ixb=ixi;  iyb=iyi; end
        ixi(ii)=[]; iyi(ii)=[]; % outside domain, note need to make sure particle not on land
        nixi=length(ixi);
        rs=zeros(nixi,1);
        w=zeros(nixi,1);
        for i=1:nixi
            ixr = ixi(i)-ix1+1+di-sum(ncutx);   % convert to local r index
            iyr = iyi(i)-iy1+1+di-sum(ncuty);
            rs(i)=r(ixr,iyr);
            %w(i) = prod(rs( i~=1:nixi) );      %Edit - RAM 1/15/21
        end
        for i=1:nixi
            w(i) = prod(rs( i~=1:nixi) );       %Edit - RAM 1/15/21
        end
        w = w/sum(w);
        v=0;      
        for i=1:nixi
            v=v+w(i)*vbar(ixi(i),iyi(i));
        end
        
        %----Now integrate with Adams-Bashforth step-----------------------
        if ~bounce
            if it == 1          %Forward step for it=1
                x = x + u*dt/1000;
                y = y + v*dt/1000;
            elseif floor(it/6)==ceil(it/6) %For times with new particles
                if ip > npold   %Forward step for new particles
                    x = x + u*dt/1000;
                    y = y + v*dt/1000;
                else            %Adams-Bashforth step
                    %x = x + (dt/(2*1000))*((3*u) - upre(ip));
                    %y = y + (dt/(2*1000))*((3*v) - vpre(ip));
                    x = x + u*dt/1000;
                    y = y + v*dt/1000;
                end
            else                %Adams-Bashforth step
                %x = x + (dt/(2*1000))*((3*u) - upre(ip));
                %y = y + (dt/(2*1000))*((3*v) - vpre(ip));
                x = x + u*dt/1000;
                y = y + v*dt/1000;
            end
        else
            if bouncu
                xb=xr(ixb(~iib),iyb(~iib));
                yb=yr(ixb(~iib),iyb(~iib));
                xix=xr(ixb(1),iyb(1));
                yix=yr(ixb(1),iyb(1));
            elseif sum(iib) == 5
                ixibig=[ix,ix+1,ix-1,ix,ix+1,ix-1,ix,ix+1,ix-1]; 
                iyibig=[iy,iy,iy,iy+1,iy+1,iy+1,iy-1,iy-1,iy-1];
                iibig = ixibig>nxu | iyibig>nyv | ixibig<1 | iyibig<1; % remove outside grid points
                ixibig(iibig)=[]; iyi(iibig)=[];
                tmp = diag(ubar(ixibig,iyibig));
                iibig = isnan(tmp);
                xb=xr(ixibig(~iibig),iyibig(~iibig));
                yb=yr(ixibig(~iibig),iyibig(~iibig));
                xix=xr(ixibig(1),iyibig(1));
                yix=yr(ixibig(1),iyibig(1));
            else
                xb=xr(ixb(~iib),iyb(~iib));
                yb=yr(ixb(~iib),iyb(~iib));
                xix=xr(ixb(1),iyb(1));
                yix=yr(ixb(1),iyb(1));
            end
            xb=diag(xb); xb = xb(1);
            dxb=xix-xb;
            x = x -0.3*dxb;
            yb=diag(yb); yb = yb(1);
            dyb=yix-yb;
            y = y -0.5*dyb;
            bounce=0;
        end
        
        %Stochastic Forcing Code-------------------------------------------
        if sto == 1 
            if it <= 8
                xs((it-savcorrect+1),ip)=x; ys((it-savcorrect+1),ip)=y;
            else
                if xs((it-savcorrect),ip) <= 180
                    xs((it-savcorrect+1),ip)=x; ys((it-savcorrect+1),ip)=y;
                else
                    walkx = randn*sqrt(edD*dt); walky = randn*sqrt(edD*dt);
                    xs((it-savcorrect+1),ip)=x + walkx; ys((it-savcorrect+1),ip)=y + walky; %With Random Walk
                end
            end
        else
            if partDrift == 1
                if sitarray(ip) > it
                    xs((it-savcorrect+1),ip)=xs((it-savcorrect),ip); ys((it-savcorrect+1),ip)=ys((it-savcorrect),ip);
                elseif eitarray(ip) < it
                    xs((it-savcorrect+1),ip)=xs((it-savcorrect),ip); ys((it-savcorrect+1),ip)=ys((it-savcorrect),ip);
                else
                    xs((it-savcorrect+1),ip)=x; ys((it-savcorrect+1),ip)=y;
                end
            else
                xs((it-savcorrect+1),ip)=x; ys((it-savcorrect+1),ip)=y;
            end
        end
        
        %Save values for later time steps for particle ip------------------
        ixs(ip)=ix; iys(ip)=iy;
        upre(ip) = u;
        vpre(ip) = v;
        
    end %ip-loop-----------------------------------------------------------
    
    %Plot Histograph to save values for each time step---------------------
    figure(fig2);
    hplot2 = histogram2(xs((it-savcorrect+1),:),ys((it-savcorrect+1),:),[180:binsize:265],[260:binsize:405],'DisplayStyle','tile','ShowEmptyBins','on');
    binval=hplot2.Values;
    if it == 1
        binvalarray = binval;
    else
        binvalarray = cat(3,binvalarray,binval);
    end
    
    %-- plot results, note isv flag for saving images for movie------------
    thsmooth=round((t(it)-t0_month)*24*10000)/10000;
    if mod(thsmooth,1.)==0
        figure(fig1);
        if hashtagbinning == 1
            h=axes('Position',[.03 .05 .45 .95], 'Layer', 'bottom');
        end
        sbar(250:300,1:60)=nan;
        pcolor(xr,yr,sbar); shading('flat')
        if plastic == 1
            colormap(h,upside_down_cmap)
            caxis([0 35])
            c1=colorbar('westoutside','YTick',[0:5:35]);
            c1.TickLabels = num2cell(0:0.5:3.5);
            colormap(c1,cmap);
            c1.Color = 0.7*[1 1 1];
            c1.Label.String = 'Particle Concentration (piece m^{-3})';
        else 
%             colormap(h,parula)
%             caxis([20 32])
            colormap(h,cmapsalt)
            caxis([0 35])
            c1=colorbar('westoutside');
            c1.Color = 0.7*[1 1 1];
            c1.Label.String = 'Salinity (PSU)';
        end
        hold on
        if partialpart == 1
            plot(xs((it-savcorrect),:),ys((it-savcorrect),:),'r.','MarkerSize',0.05)
        elseif partDrift == 1
            plot(xs((it-savcorrect),:),ys((it-savcorrect),:),'k.','MarkerSize',5)
        else
            plot(xs((it-savcorrect),:),ys((it-savcorrect),:),'k.','MarkerSize',0.05)
        end
        set(gca,'DataAspectRatio',[1 1 1]) 
        set(gca,'Color','k')
        set(gca,'XColor',0.7*[1 1 1],'YColor',0.7*[1 1 1])
        set(gca,'XLim',[180  260],'YLim', [260  400])   %x=[180 260] y=[260  400]
        %set(gca,'XLim',[170  250],'YLim', [260  400])   %x=[180 260] y=[260  400]
        set(gca,'FontSize',12)
        xlabel('x (km)')
        xticks([180,200,220,240,260])
        %xticks([170,190,210,230,250])
        yticks([260,280,300,320,340,360,380,400])
        ylabel('y (km)')
        if hashtagbinning == 0
            text(210,390,datestr(t(it)),'fontsize',16,'Color',0.7*[1 1 1])
        end
        
        %%%Code for Binning Data-------------------------------------------
        if hashtagbinning == 1
            %%%bin size visual
            if binshow == 1
                if it <= 12     %Multiples of 6
                    for bininy = 290:binsize:400
                        refline(0,bininy);
                    end
                    for bininx = 180:binsize:245
                        line([bininx bininx],[260 400]);
                    end
                end
            end
            %title({['Number of Particles: ',num2str(np)]},'Color',0.7*[1 1 1])
            %text(210,390,datestr(t(it)),'fontsize',14.5,'Color',0.0*[1 1 1]) %%%Salinity Plot title
                        
            %xsdos=(xs((it-savcorrect+1),:))/np; ysdos=(ys((it-savcorrect+1),:))/np;
            
            h2=axes('Position',[.57 .11 .38 .82],'Layer', 'top');
            set(gca,'XColor',0.7*[1 1 1],'YColor',0.7*[1 1 1],'ZColor',0.7*[1 1 1])
            set(gca, 'zscale','log');
            %hplot2 = histogram2(xs((it-savcorrect+1),:),ys((it-savcorrect+1),:),[180:binsize:265],[260:binsize:405],'DisplayStyle','tile','ShowEmptyBins','on');
            %binvalcheck=hplot2.Values;
            %set(gca, 'zscale','log');
            if usenormbin == 1
                binadj = binval/binvalmax;
                pcolor(y6,x6,binadj); shading('flat');
                %caxis([0.000001 .1])
                caxis([0.00001 1])
            else
                hplot2 = histogram2(xs((it-savcorrect+1),:),ys((it-savcorrect+1),:),[180:binsize:265],[260:binsize:405],'DisplayStyle','tile','ShowEmptyBins','on');
                %hplot2 = histogram2(xs((it-savcorrect+1),:),ys((it-savcorrect+1),:),[170:binsize:255],[260:binsize:405],'DisplayStyle','tile','ShowEmptyBins','on');
                caxis([0.1 10000])
            end
            colormap(h2,mymap)
            cb=colorbar;
            set(h2,'ColorScale','log')
            %caxis([0 .001])
            %caxis([0.000001 .1])
            %caxis([0.00001 1])
            cb.Color = 0.7*[1 1 1];
            cb.Label.String = 'Number of Particles';
            %caxis([1 300]);
            set(gca,'XColor',0.7*[1 1 1],'YColor',0.7*[1 1 1],'ZColor',0.7*[1 1 1])
            xlim([180 260])
            ylim([260 400])
            %set(gca,'XLim',[170  250],'YLim', [260  400])
            xlabel('x (km)')
            ylabel('y (km)')
            zlabel('Number of Particles')
            set(gca,'FontSize',12)
            %title(['Size of Bin in square km: ',num2str((binsize)^2)],'Color',0.7*[1 1 1])
            hold on
            plot(ROMSshore,'FaceColor','k','FaceAlpha',.95)
%             plot(lowpolyx,lowpolyy,'LineWidth',2,'Color','c');
%             plot(uppolyx,uppolyy,'LineWidth',2,'Color','c');
            %text(210,390,datestr(t(it)),'fontsize',14.5,'Color',0.7*[1 1 1]) %%%Salinity Plot title
            
            timestamp = itime-0;    %Correction Factors: -0=UTC -4=EDT -5=EST
            if itime == 730
                otc = (timeold/3600/24)+t0_ROMS;
            elseif timestamp < 1
                timestamp = timestamp+730;
                otc = (timeold/3600/24)+t0_ROMS;
            else
                otc = (time/3600/24)+t0_ROMS;
            end
            text(202,390,datestr(otc(timestamp)),'fontsize',14.5,'Color',0.7*[1 1 1]) %%%Salinity Plot title
            %text(198,390,[datestr(otc(timestamp)),' EDT'],'fontsize',14.5,'Color',0.7*[1 1 1]) %%%Salinity Plot title
            %text(189,393,[datestr(otc(timestamp)),' EDT'],'fontsize',14.5,'Color',0.7*[1 1 1]) %%%Salinity Plot title
        end
        %%%----------------------------------------------------------------
      
        if isv==1
            cnt=cnt+1;
            print([dirname,num2str(cnt)],'-dpng');
            %figname = [dirname,num2str(cnt)];
            %export_fig figname -png -nocrop
        end
        
        pause(0.01)
        clf(fig1)
        
    end %End plotting if-statement
    
    if smallsave == 1 && lensave == (it-savcorrect+1)
        savnum = savnum+1;              % Update966s the number on the save file
        savname = [smallname,num2str(savnum)];  % Name for save file
        enddatetime = datestr((t(it))); % The last time recorded in the file
        save(savname,'binvalarray','enddatetime','np','nporiginal','startdatetime','xs','ys'); % Saves essential variables to be rewritten
        xslast = xs(end,:);             % Saves last value for use in the new array
        yslast = ys(end,:);
        binlast = squeeze(binvalarray(:,:,end));
        clear xs ys binvalarray         % Clears old values to clear up memory
        xs = zeros(lensave,np);         % Creates new xs array for the next memory file
        ys = zeros(lensave,np);         % Same as x
        xs(1,:) = xslast;               % Takes the last x position in the previous file and makes it the first in the new file
        ys(1,:) = yslast;               % Same as x
        binvalarray = binlast;          % Idea as x/y position, except takes the full binval array in x and y directions
        savcorrect = it;                % Updates the savcorrect counter so that way (it-savcorrect) counts from 1-72
        startdatetime = enddatetime;    % Due to overlap of one hour, the first time of the new file is the same as the last time of the previous
    end
    
    %pause(0.01)
    
end %End time loop

%%
%%%Saves last xs/ys/binvalarray files--------------------------------------
savnum = savnum+1;                      % Updates the number on the save file
savname = [smallname,num2str(savnum)];  % Name for save file
enddatetime = datestr((t(it)));         % The last time recorded in the file
xs(all(~xs,2),:)=[];                    % Cuts out any extra rows (rows of zeros)
ys(all(~ys,2),:) = [];                  % Same as x
save(savname,'binvalarray','enddatetime','np','nporiginal','startdatetime','xs','ys'); % Saves essential variables to be rewritten
%--------------------------------------------------------------------------

toc

