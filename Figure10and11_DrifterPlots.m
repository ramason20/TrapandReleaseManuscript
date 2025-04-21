%%%Plot Drifter Tracks based on Reformatted, Generic Name Drifters put on
%%%the drifter repository in the Google Drive
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 12/8/23 
%Edited: 12/8/23
%------------------------------------------------------------------------

%%%LAST EDIT - Change to Generic Drifter Names

tic

%% - Flags---------------------------------------------------------------
isv = 0;
plotmean = 0;

%% - Load Data-----------------------------------------------------------

%June22Release - Use drifters 1-10
    load('\savloc\2022Jun14Release_TimeInterp.mat');
    t0EDT = datetime(2022,06,14,11,00,0);
    tendEDT = datetime(2022,06,15,11,00,0);
    tdrift = t0EDT:minutes(5):tendEDT;
    endpoint = 169;
    drifterloncat = cat(2,drift2_lon,drift3_lon,drift4_lon,drift5_lon,drift6_lon,drift7_lon,drift8_lon,drift9_lon);
    drifterlatcat = cat(2,drift2_lat,drift3_lat,drift4_lat,drift5_lat,drift6_lat,drift7_lat,drift8_lat,drift9_lat);
    %drifterloncat = cat(2,drift2_lon,drift3_lon,drift4_lon,drift5_lon,drift6_lon);
    %drifterlatcat = cat(2,drift2_lat,drift3_lat,drift4_lat,drift5_lat,drift6_lat);
    Jun22_driftmeanlon = mean(drifterloncat,2);     Jun22_driftmeanlon = squeeze(Jun22_driftmeanlon(1:endpoint,:));
    Jun22_driftmeanlat = mean(drifterlatcat,2);     Jun22_driftmeanlat = squeeze(Jun22_driftmeanlat(1:endpoint,:));

%July23Release - Use drifters 1-6
    load('\savloc\2023Jul27Release_TimeInterp.mat');
    t0EDT = datetime(2023,07,27,9,00,0);
    tendEDT = datetime(2023,07,31,11,00,0);
    tdrift = t0EDT:minutes(5):tendEDT;
    endpoint = 163;
    %drifterloncat = cat(2,drift1_lon,drift2_lon,drift3_lon,drift4_lon,drift5_lon,drift6_lon,drift7_lon,drift8_lon,drift9_lon,drift10_lon,drift11_lon,drift12_lon);
    %drifterlatcat = cat(2,drift1_lat,drift2_lat,drift3_lat,drift4_lat,drift5_lat,drift6_lat,drift7_lat,drift8_lat,drift9_lat,drift10_lat,drift11_lat,drift12_lat);
    drifterloncat = cat(2,drift1_lon,drift2_lon,drift3_lon,drift5_lon,drift6_lon);
    drifterlatcat = cat(2,drift1_lat,drift2_lat,drift3_lat,drift5_lat,drift6_lat);
    Jul23_driftmeanlon = mean(drifterloncat,2);     Jul23_driftmeanlon = squeeze(Jul23_driftmeanlon(1:endpoint,:));
    Jul23_driftmeanlat = mean(drifterlatcat,2);     Jul23_driftmeanlat = squeeze(Jul23_driftmeanlat(1:endpoint,:));
    
%Aug23Release - Use drifters 1-6
    load('\savloc\2023Aug01Release_TimeInterp.mat');
    dirname = 'F:\RAM_Research\Drifter_Data\2023Aug\AugDrifter_Vid\';
    t0EDT = datetime(2023,08,01,10,00,0);
    tendEDT = datetime(2023,08,02,11,00,0);
    tdrift = t0EDT:minutes(5):tendEDT;
    endpoint = 151;
    drifterloncat = cat(2,drift1_lon,drift2_lon,drift3_lon,drift4_lon,drift5_lon,drift6_lon);
    drifterlatcat = cat(2,drift1_lat,drift2_lat,drift3_lat,drift4_lat,drift5_lat,drift6_lat);
    Aug23_driftmeanlon = mean(drifterloncat,2);     Aug23_driftmeanlon = squeeze(Aug23_driftmeanlon(1:endpoint,:));
    Aug23_driftmeanlat = mean(drifterlatcat,2);     Aug23_driftmeanlat = squeeze(Aug23_driftmeanlat(1:endpoint,:));

load('\savloc\NOAA_HighResBathy.mat');
cmap = cmocean('ice',25);

%% - Generate Figure-----------------------------------------------------
fig1=figure;                                    % Used for plotting main figure
set(gcf,'PaperSize',[13.33/2 7.5])              % Sets scale of figure so it does not become distorted
%set(gcf, 'InvertHardcopy', 'off')               % Tells Matlab to save the figure background color and not to change it to white
set(gcf,'Position',[378   302   600   600])     % Sets location and size of the plot
cnt = 1000;                                     % Imagine number counting for saving figures

% - Plot Loop-----------------------------------------------------------

%for itime = iloop
    figure(fig1);
    pcolor(loni,lati,HRbathy);  shading('flat');
    hold on

    %Build Track
    if plotmean == 1
        plot(Jul23_driftmeanlon,Jul23_driftmeanlat,'r','LineWidth',2);
        plot(Aug23_driftmeanlon,Aug23_driftmeanlat,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
        plot(Jun22_driftmeanlon,Jun22_driftmeanlat,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
    else
        plot(drift1_lon(1:endpoint),drift1_lat(1:endpoint),'w'); 
        plot(drift2_lon(1:endpoint),drift2_lat(1:endpoint),'w');
        plot(drift3_lon(1:endpoint),drift3_lat(1:endpoint),'w');
        plot(drift4_lon(1:endpoint),drift4_lat(1:endpoint),'w'); 
        plot(drift5_lon(1:endpoint),drift5_lat(1:endpoint),'w');
        plot(drift6_lon(1:endpoint),drift6_lat(1:endpoint),'w');
        plot(drift7_lon(1:endpoint),drift7_lat(1:endpoint),'w');
        plot(drift8_lon(1:endpoint),drift8_lat(1:endpoint),'w');
        plot(drift9_lon(1:endpoint),drift9_lat(1:endpoint),'w');
        plot(drift10_lon(1:endpoint),drift10_lat(1:endpoint),'w'); 
      % plot(Jul23_driftmeanlon,Jul23_driftmeanlat,'r','LineWidth',2);
      % plot(Aug23_driftmeanlon,Aug23_driftmeanlat,'r','LineWidth',2);
        plot(Jun22_driftmeanlon,Jun22_driftmeanlat,'r','LineWidth',2);
    end
    
    %Initial Position and Current Position
    if plotmean == 1
        plot(Jul23_driftmeanlon(1),Jul23_driftmeanlat(1),'rs','MarkerFaceColor','w');      plot(Jul23_driftmeanlon(end),Jul23_driftmeanlat(end),'rp','MarkerFaceColor','k','MarkerSize',8);
        plot(Aug23_driftmeanlon(1),Aug23_driftmeanlat(1),'Color',[0.9290 0.6940 0.1250],'Marker','square','MarkerFaceColor','w');      plot(Aug23_driftmeanlon(end),Aug23_driftmeanlat(end),'Color',[0.9290 0.6940 0.1250],'Marker','p','MarkerFaceColor','k','MarkerSize',8);
        plot(Jun22_driftmeanlon(1),Jun22_driftmeanlat(1),'Color',[0.9290 0.6940 0.1250],'Marker','square','MarkerFaceColor','w');      plot(Jun22_driftmeanlon(end),Jun22_driftmeanlat(end),'Color',[0.9290 0.6940 0.1250],'Marker','p','MarkerFaceColor','k','MarkerSize',8);
    else
        plot(drift1_lon(1),drift1_lat(1),'ws','MarkerFaceColor','w');      plot(drift1_lon(endpoint),drift1_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift2_lon(1),drift2_lat(1),'ws','MarkerFaceColor','w');      plot(drift2_lon(endpoint),drift2_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift3_lon(1),drift3_lat(1),'ws','MarkerFaceColor','w');      plot(drift3_lon(endpoint),drift3_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift4_lon(1),drift4_lat(1),'ws','MarkerFaceColor','w');      plot(drift4_lon(endpoint),drift4_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift5_lon(1),drift5_lat(1),'ws','MarkerFaceColor','w');      plot(drift5_lon(endpoint),drift5_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift6_lon(1),drift6_lat(1),'ws','MarkerFaceColor','w');      plot(drift6_lon(endpoint),drift6_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift7_lon(1),drift7_lat(1),'ws','MarkerFaceColor','w');      plot(drift7_lon(endpoint),drift7_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift8_lon(1),drift8_lat(1),'ws','MarkerFaceColor','w');      plot(drift8_lon(endpoint),drift8_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift9_lon(1),drift9_lat(1),'ws','MarkerFaceColor','w');      plot(drift9_lon(endpoint),drift9_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
        plot(drift10_lon(1),drift10_lat(1),'ws','MarkerFaceColor','w');    plot(drift10_lon(endpoint),drift10_lat(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
%         plot(Jul23_driftmeanlon(1),Jul23_driftmeanlat(1),'rs','MarkerFaceColor','w');      plot(Jul23_driftmeanlon(end),Jul23_driftmeanlat(end),'rp','MarkerFaceColor','k','MarkerSize',8);
%         plot(Aug23_driftmeanlon(1),Aug23_driftmeanlat(1),'rs','MarkerFaceColor','w');      plot(Aug23_driftmeanlon(end),Aug23_driftmeanlat(end),'rp','MarkerFaceColor','k','MarkerSize',8);
        plot(Jun22_driftmeanlon(1),Jun22_driftmeanlat(1),'rs','MarkerFaceColor','w');      plot(Jun22_driftmeanlon(end),Jun22_driftmeanlat(end),'rp','MarkerFaceColor','k','MarkerSize',8);
    end

    colormap(cmap)
    c1 = colorbar;
    caxis([-20 0])
    c1.Label.String = 'Depth (m)';                      %Colorbar label
    
    set(gca,'DataAspectRatio',[1 1 1])                  %Adjusts aspect ration of figure
    %set(gca,'XLim',[-75.5  -75.1],'YLim', [38.9  39.4])  %Full Bay
    set(gca,'XLim',[-75.35  -75.1],'YLim', [38.925  39.175])  %Lower Bay
            plot(drift1_lon(61),drift1_lat(61),'yo','MarkerFaceColor','y'); 
            plot(drift2_lon(61),drift2_lat(61),'yo','MarkerFaceColor','y');
            plot(drift3_lon(61),drift3_lat(61),'yo','MarkerFaceColor','y');
            plot(drift4_lon(61),drift4_lat(61),'yo','MarkerFaceColor','y'); 
            plot(drift5_lon(61),drift5_lat(61),'yo','MarkerFaceColor','y');
            plot(drift6_lon(61),drift6_lat(61),'yo','MarkerFaceColor','y');
            plot(drift7_lon(61),drift7_lat(61),'yo','MarkerFaceColor','y'); 
            plot(drift8_lon(61),drift8_lat(61),'yo','MarkerFaceColor','y');
            plot(drift9_lon(61),drift9_lat(61),'yo','MarkerFaceColor','y');
            plot(drift10_lon(61),drift10_lat(61),'yo','MarkerFaceColor','y');
            load('\savloc\CTDData_2022Jun14.mat');
            lattran = tran3_lat;    lontran = -tran3_lon;
%     set(gca,'XLim',[-75.5  -75.3],'YLim', [39.2  39.4])  %Upper Bay - Spring
%             plot(drift1_lon(57),drift1_lat(57),'yo','MarkerFaceColor','y'); 
%             plot(drift2_lon(57),drift2_lat(57),'yo','MarkerFaceColor','y');
%             plot(drift3_lon(57),drift3_lat(57),'yo','MarkerFaceColor','y');
%             plot(drift4_lon(57),drift4_lat(57),'yo','MarkerFaceColor','y'); 
%             plot(drift5_lon(57),drift5_lat(57),'yo','MarkerFaceColor','y');
%             plot(drift6_lon(57),drift6_lat(57),'yo','MarkerFaceColor','y');
%             load('\savloc\UpperBaySpringTideCTD.mat');
%             tranno = 2;     transectfull = 2:9;
%             tranf = transectfull(1);    trane = transectfull(end);
%             lattran = squeeze(lat(tranf:trane));    lontran = squeeze(lon(tranf:trane));
%     set(gca,'XLim',[-75.5  -75.325],'YLim', [39.275  39.4])  %Upper Bay - Neap
%             plot(drift1_lon(121),drift1_lat(121),'yo','MarkerFaceColor','y'); 
%             plot(drift2_lon(121),drift2_lat(121),'yo','MarkerFaceColor','y');
%             plot(drift3_lon(121),drift3_lat(121),'yo','MarkerFaceColor','y');
%             plot(drift5_lon(121),drift5_lat(121),'yo','MarkerFaceColor','y');
%             plot(drift6_lon(121),drift6_lat(121),'yo','MarkerFaceColor','y');
%             load('\savloc\UpperBayNeapTideCTD.mat');
%             tranno = 16;    transectfull = 770:823;
%             tranf = transectfull(1);    trane = transectfull(end);
%             lattran = squeeze(lat(tranf:trane));    lontran = squeeze(lon(tranf:trane));
    xticks([-75.5,-75.45,-75.4,-75.35,-75.3,-75.25,-75.2,-75.15,-75.1,-75.05,-75.0,-74.95,-74.9])
    yticks([38.9,38.95,39.0,39.05,39.1,39.15,39.2,39.25,39.3,39.35,39.4,39.45,39.5,39.55,39.6])
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'FontSize',12)
    title('Spring Tide Drifter Tracks')
    %title('Neap Tide Drifter Tracks')
    grid on

    plot(lontran,lattran,'y-','LineWidth',2);

     if isv == 1
         if plotmean == 1
             print('\savloc\TrapReleaseFigure_MeanPlot','-dpng')
         else
             print('\savloc\TrapReleaseFigure_SpringNew','-dpng')
         end
     end

%% - Zoom Plot

%Convert from xr/yr into lon/lat
xrs = reshape(xbathy,[3397*2088, 1]);
yrs = reshape(ybathy,[3397*2088, 1]);
lonrs = reshape(loni,[3397*2088, 1]);
latrs = reshape(lati,[3397*2088, 1]);
Fx = scatteredInterpolant(lonrs,latrs,xrs);
Fy = scatteredInterpolant(lonrs,latrs,yrs);

drift1x = Fx(drift1_lon(1:endpoint),drift1_lat(1:endpoint));     drift1y = Fy(drift1_lon(1:endpoint),drift1_lat(1:endpoint));
drift2x = Fx(drift2_lon(1:endpoint),drift2_lat(1:endpoint));     drift2y = Fy(drift2_lon(1:endpoint),drift2_lat(1:endpoint));
drift3x = Fx(drift3_lon(1:endpoint),drift3_lat(1:endpoint));     drift3y = Fy(drift3_lon(1:endpoint),drift3_lat(1:endpoint));
drift4x = Fx(drift4_lon(1:endpoint),drift4_lat(1:endpoint));     drift4y = Fy(drift4_lon(1:endpoint),drift4_lat(1:endpoint));
drift5x = Fx(drift5_lon(1:endpoint),drift5_lat(1:endpoint));     drift5y = Fy(drift5_lon(1:endpoint),drift5_lat(1:endpoint));
drift6x = Fx(drift6_lon(1:endpoint),drift6_lat(1:endpoint));     drift6y = Fy(drift6_lon(1:endpoint),drift6_lat(1:endpoint));
drift7x = Fx(drift7_lon(1:endpoint),drift7_lat(1:endpoint));     drift7y = Fy(drift7_lon(1:endpoint),drift7_lat(1:endpoint));
drift8x = Fx(drift8_lon(1:endpoint),drift8_lat(1:endpoint));     drift8y = Fy(drift8_lon(1:endpoint),drift8_lat(1:endpoint));
drift9x = Fx(drift9_lon(1:endpoint),drift9_lat(1:endpoint));     drift9y = Fy(drift9_lon(1:endpoint),drift9_lat(1:endpoint));
drift10x = Fx(drift10_lon(1:endpoint),drift10_lat(1:endpoint));  drift10y = Fy(drift10_lon(1:endpoint),drift10_lat(1:endpoint));
%Jul23x = Fx(Jul23_driftmeanlon,Jul23_driftmeanlat);              Jul23y = Fy(Jul23_driftmeanlon,Jul23_driftmeanlat);
%Aug23x = Fx(Aug23_driftmeanlon,Aug23_driftmeanlat);              Aug23y = Fy(Aug23_driftmeanlon,Aug23_driftmeanlat);
Jun22x = Fx(Jun22_driftmeanlon,Jun22_driftmeanlat);              Jun22y = Fy(Jun22_driftmeanlon,Jun22_driftmeanlat);


%% - Generate Figure-----------------------------------------------------
fig1=figure;                                    % Used for plotting main figure
set(gcf,'PaperSize',[13.33/2 7.5])              % Sets scale of figure so it does not become distorted
%set(gcf, 'InvertHardcopy', 'off')               % Tells Matlab to save the figure background color and not to change it to white
set(gcf,'Position',[378   302   600   600])     % Sets location and size of the plot
cnt = 1000;                                     % Imagine number counting for saving figures

% - Plot Loop-----------------------------------------------------------

%for itime = iloop
    figure(fig1);
    pcolor(xbathy,ybathy,HRbathy);  shading('flat');
    hold on

    plot(drift1x,drift1y,'w'); 
    plot(drift2x,drift2y,'w');
    plot(drift3x,drift3y,'w');
    plot(drift4x,drift4y,'w');
    plot(drift5x,drift5y,'w');
    plot(drift6x,drift6y,'w');
%     plot(drift7x,drift7y,'w');
%     plot(drift8x,drift8y,'w');
%     plot(drift9x,drift9y,'w');
%     plot(drift10x,drift10y,'w');
    plot(Jul23x,Jul23y,'r','LineWidth',2);
    %plot(Aug23x,Aug23y,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
%     plot(Jun22x,Jun22y,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);
    
    %Initial Position and Current Position
    plot(drift1x(1),drift1y(1),'ws','MarkerFaceColor','w');      plot(drift1x(endpoint),drift1y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
    plot(drift2x(1),drift2y(1),'ws','MarkerFaceColor','w');      plot(drift2x(endpoint),drift2y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
    plot(drift3x(1),drift3y(1),'ws','MarkerFaceColor','w');      plot(drift3x(endpoint),drift3y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
    plot(drift4x(1),drift4y(1),'ws','MarkerFaceColor','w');      plot(drift4x(endpoint),drift4y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
    plot(drift5x(1),drift5y(1),'ws','MarkerFaceColor','w');      plot(drift5x(endpoint),drift5y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
    plot(drift6x(1),drift6y(1),'ws','MarkerFaceColor','w');      plot(drift6x(endpoint),drift6y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
%     plot(drift7x(1),drift7y(1),'ws','MarkerFaceColor','w');      plot(drift7x(endpoint),drift7y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
%     plot(drift8x(1),drift8y(1),'ws','MarkerFaceColor','w');      plot(drift8x(endpoint),drift8y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
%     plot(drift9x(1),drift9y(1),'ws','MarkerFaceColor','w');      plot(drift9x(endpoint),drift9y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
%     plot(drift10x(1),drift10y(1),'ws','MarkerFaceColor','w');    plot(drift10x(endpoint),drift10y(endpoint),'kp','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
%     plot(Jul23x(1),Jul23y(1),'rs','MarkerFaceColor','w');        plot(Jul23x(end),Jul23y(end),'rp','MarkerFaceColor','k','MarkerSize',8);
    plot(Aug23x(1),Aug23y(1),'Color',[0.9290 0.6940 0.1250],'Marker','square','MarkerFaceColor','w');      plot(Aug23x(end),Aug23y(end),'Color',[0.9290 0.6940 0.1250],'Marker','p','MarkerFaceColor','k','MarkerSize',8);
%     plot(Jun22x(1),Jun22y(1),'Color',[0.9290 0.6940 0.1250],'Marker','square','MarkerFaceColor','w');      plot(Jun22x(end),Jun22y(end),'Color',[0.9290 0.6940 0.1250],'Marker','p','MarkerFaceColor','k','MarkerSize',8);

    colormap(cmap)
    c1 = colorbar;
    caxis([-20 0])
    c1.Label.String = 'Depth (m)';                      %Colorbar label
    
    set(gca,'DataAspectRatio',[1 1 1])                  %Adjusts aspect ration of figure
    xlabel('x (km)')
    ylabel('y (km)')
    set(gca,'FontSize',12)
    title('Lower Bay Spring Tide Drifter Tracks')
    %xlim([175 195]);     ylim([367 392])  %Upper Spring
    xlim([190 215]);     ylim([325 360])  %Lower Spring
    %title('Neap Tide Drifter Tracks')
    %xlim([174 191]);     ylim([376 392])
    grid on

    if isv == 1
        print('\savloc\TrapReleaseFigure_LowerSpring_Zoom','-dpng')
    end

%toc