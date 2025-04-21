%%%Plots Residence Time line with 3 Forcings for Idealized Cases
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 4/2/20 
%Edited: 4/21/25
%------------------------------------------------------------------------

tic

%%%Flag to select which subset of particles to use
AllPart = 1;
addpart = 0;

%%%Sets up parameters
dirname = '\ParticleLocation\PartTrack\';
savfigname = '\savloc\';
MarchRST = 1;
Jan08RST = 0;
numofiles = 21; %42; %15; %42; %57;  %Total Number of PartTrack Files
killline=0;      %Use Kill Line/Do not allow particles to reenter the bay?
binsize = 10;    %Sets size of histogram bins
mymap = [0/255 0/255 128/255; 0/255 128/255 255/255; 46/255 151/255 208/255; 92/255 174/255 162/255; 139/255 197/255 115/255; ...
        185/255 220/255 69/255; 231/255 243/255 23/255; 255/255 243/255 0/255; 255/255 220/255 0/255; 255/255 197/255 0/255; ...
        255/255 174/255 0/255; 255/255 151/255 0/255; 255/255 128/255 0/255; 255/255 96/255 0/255; 255/255 64/255 0/255; ...
        255/255 32/255 0/255; 255/255 0/255 0/255; 170/255 85/255 0/255; 85/255 170/255 0/255; 0/255 255/255 0/255];

%Creation of Kill Line
xDE = 220;  yDE = 298;  %First values used: x=215 y=295
xNJ = 250;  yNJ = 334;  %First values used: x=245 y=330
numopoints = 1000;
outxline=linspace(xDE,xNJ,numopoints);
outyline=linspace(yDE,yNJ,numopoints);

%Loads in parameters
load('\ParticleStartUpFile.mat')
clear excelpart

rowstart = 145:143:1432; rowstart = cat(2,1,rowstart);
rowend = 287:143:1432;   rowend = cat(2,144,rowend,1440);

for ifile = 1:numofiles
    %%%Load in data
    load([dirname,num2str(ifile),'.mat']);
    clear binvalarray
    if AllPart == 1
        if addpart == 1
            if ifile >1
                xs = squeeze(xs(2:end,:));
                ys = squeeze(ys(2:end,:));
            end
            npseg = 55614; %50534;
            xsbig = ones(size(xs,1),npseg).*190;
            ysbig = ones(size(ys,1),npseg).*380;
            xsbig(:,1:size(xs,2)) = xs;
            ysbig(:,1:size(xs,2)) = ys;
            xs = xsbig; ys = ysbig;
        else
            if ifile > 1
                xs=squeeze(xs(2:end,:));
                ys=squeeze(ys(2:end,:));
            end
            npseg=size(xs,2);
        end
    else
        xs=squeeze(xs(:,xseg));
        ys=squeeze(ys(:,yseg));
        npseg=length(xseg);
    end
    begtime=datenum(startdatetime);
    
    %%%Sets up needed arrays
    lentime = size(xs,1); 
    timeleft = zeros(1,npseg);
    notinbay = zeros(lentime,npseg);

    %%%Sets up needed parameters
    sdt = datenum(startdatetime);   %Units are days
    finalspotx = zeros(1,npseg);
    finalspoty = zeros(1,npseg);
    xbin = 180:binsize:260;
    ybin = 260:binsize:400;
    
    %%%Determines if particle left bay
    for it = 1:lentime

        modtime = sdt+(it/24/12);   %Changes it from increments of 5 minutes to days

        for ip = 1:npseg
            txid = xs(it,ip);
            tyid = ys(it,ip); 
            for index = 2:numopoints
                xbo=outxline(index);
                ybo=outyline(index-1);
                if txid >= 250
                    notinbay(it,ip)=1;
                elseif tyid <= 290
                    notinbay(it,ip)=1;
                elseif txid >= xbo
                    if tyid <= ybo
                        notinbay(it,ip)=1;
                    end
                end
            end
            
        end
    end

    if ifile == 1
        notinbayfull = notinbay;
    else
        notinbayfull = cat(1,notinbayfull,notinbay);
    end
    
end

fulltime = size(notinbayfull,1);
leftbay = zeros(fulltime,npseg);

for itleft = 1:fulltime
    for ip2 = 1:npseg
        if itleft >= 4
            if killline == 1
                if leftbay((itleft-1),ip2) == 1
                    leftbay(itleft,ip2) = 1;
                elseif notinbayfull(itleft,ip2) == 0
                    %Nothing to be done
                elseif notinbayfull((itleft-3),ip2)==1 && notinbayfull((itleft-2),ip2)==1 && notinbayfull((itleft-1),ip2)==1 && notinbayfull(itleft,ip2)==1
                    leftbay(itleft,ip2) = 1;
                end
            else
                if notinbayfull(itleft,ip2) == 0
                    %Nothing to be done
                elseif notinbayfull((itleft-3),ip2)==1 && notinbayfull((itleft-2),ip2)==1 && notinbayfull((itleft-1),ip2)==1 && notinbayfull(itleft,ip2)==1
                    leftbay(itleft,ip2) = 1;
                end
            end
        end
    end
end
toc
        
%% - Time Series of when particles leave the bay
        
    numpart=npseg*ones(1,fulltime);
    xtrange = 1:(12*24):fulltime;
    for it3 = 1:fulltime
        for ip3 = 1:npseg
            inbayq = leftbay(it3,ip3);
            if inbayq == 1
                numpart(it3) = numpart(it3)-1;
            end
        end
    end
    numfrac=numpart./npseg;
    numper=numfrac.*100;
          
    %Plot Code
    fig2=figure;
    set(gcf,'Position',[250   250   1400   400])
    plot(numfrac,'k')
    xlabel('Day')
    ylabel('Fraction of Particles');
    xlim([1 fulltime]);
    xticks(xtrange) %5 minute step
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    ylim([0 1])
    yticks(0:0.1:1)
    title('Fraction of Particles in the Bay - Spring Tide')
    grid on
    set(gca,'FontSize',16)
    
    print([savfigname,'ParticlePercentage_New'],'-dpng');

toc