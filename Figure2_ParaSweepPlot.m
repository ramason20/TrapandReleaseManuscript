%%% Figure 2 Plot Code for Surface Convergence with Varying Vertical Eddy
%%% Viscosity
%------------------------------------------------------------------------
%Author: R Alan Mason
%Created: 9/19/24 
%Edited: 4/21/25
%------------------------------------------------------------------------

% - Last Edit = Initial Creation

coriolis = 1;
isv = 1;

%% - Load Data

if coriolis == 1

    rootdir = '\dataloc\';
    
    %Akv = 0.5e-3
    load([rootdir,'Akv_05\ModelNumerics\rk_suvw_nx151B10000h12Km50g099994sx200it20000000.mat']);
    akv05 = u;      akv05 = squeeze(akv05(:,76));
    akv05v = v;     akv05v = squeeze(akv05v(:,76));
    
    %Akv = 0.75e-3
    load([rootdir,'Akv_075\ModelNumerics\rk_suvw_nx151B10000h12Km75g099994sx200it20000000.mat']);
    akv075 = u;     akv075 = squeeze(akv075(:,76));
    akv075v = v;    akv075v = squeeze(akv075v(:,76));
    
    %Akv = 1.0e-3
    load([rootdir,'Akv_1\ModelNumerics\rk_suvw_nx151B10000h12Km100g099994sx200it20000000.mat']);
    akv10 = u;      akv10 = squeeze(akv10(:,76));
    akv10v = v;     akv10v = squeeze(akv10v(:,76));
    
    %Akv = 1.25e-3
    load([rootdir,'Akv_125\ModelNumerics\rk_suvw_nx151B10000h12Km125g10001sx200it20000000.mat']);
    akv125 = u;     akv125 = squeeze(akv125(:,76));
    akv125v = v;    akv125v = squeeze(akv125v(:,76));
    
    %Akv = 1.5e-3
    load([rootdir,'Akv_15\ModelNumerics\rk_suvw_nx151B10000h12Km150g1sx200it20000000.mat'])
    akv15 = u;      akv15 = squeeze(akv15(:,76));
    akv15v = v;     akv15v = squeeze(akv15v(:,76));
    
    %Akv = 1.75e-3
    load([rootdir,'Akv_175\ModelNumerics\rk_suvw_nx151B10000h12Km175g1sx200it20000000.mat']);
    akv17 = u;      akv17 = squeeze(akv17(:,76));
    akv17v = v;     akv17v = squeeze(akv17v(:,76));
    
    %Akv = 2.0e-3
    load([rootdir,'Akv_2\ModelNumerics\rk_suvw_nx151B10000h12Km200g1sx200it20000000.mat'])
    akv20 = u;      akv20 = squeeze(akv20(:,76));
    akv20v = v;     akv20v = squeeze(akv20v(:,76));

    %Akv = 2.25e-3
    load([rootdir,'Akv_225\ModelNumerics\rk_suvw_nx151B10000h12Km225g15625sx200it20000000.mat']);
    akv225 = u;     akv225 = squeeze(akv225(:,76));
    akv225v = v;    akv225v = squeeze(akv225v(:,76));
    
    %Akv = 2.5e-3
    load([rootdir,'Akv_25\ModelNumerics\rk_suvw_nx151B10000h12Km250g17361sx200it20000000.mat'])
    akv25 = u;      akv25 = squeeze(akv25(:,76));
    akv25v = v;     akv25v = squeeze(akv25v(:,76));
    
    %Akv = 2.75e-3
    load([rootdir,'Akv_275\ModelNumerics\rk_suvw_nx151B10000h12Km275g19097sx200it20000000.mat']);
    akv27 = u;      akv27 = squeeze(akv27(:,76));
    akv27v = v;     akv27v = squeeze(akv27v(:,76));
    
    %Akv = 3.0e-3
    load([rootdir,'Akv_3\ModelNumerics\rk_suvw_nx151B10000h12Km300g20833sx200it20000000.mat'])
    akv30 = u;      akv30 = squeeze(akv30(:,76));
    akv30v = v;     akv30v = squeeze(akv30v(:,76));

    %Akv = 3.25e-3
    load([rootdir,'Akv_325\ModelNumerics\rk_suvw_nx151B10000h12Km325g22569sx200it20000000.mat']);
    akv325 = u;     akv325 = squeeze(akv325(:,76));
    akv325v = v;    akv325v = squeeze(akv325v(:,76));
    
    %Akv = 3.5e-3
    load([rootdir,'Akv_35\ModelNumerics\rk_suvw_nx151B10000h12Km350g24306sx200it20000000.mat'])
    akv35 = u;      akv35 = squeeze(akv35(:,76));
    akv35v = v;     akv35v = squeeze(akv35v(:,76));
    
    %Akv = 3.75e-3
    load([rootdir,'Akv_375\ModelNumerics\rk_suvw_nx151B10000h12Km375g26042sx200it20000000.mat']);
    akv37 = u;      akv37 = squeeze(akv37(:,76));
    akv37v = v;     akv37v = squeeze(akv37v(:,76));
    
    %Akv = 4.0e-3
    load([rootdir,'Akv_4\ModelNumerics\rk_suvw_nx151B10000h12Km400g27778sx200it20000000.mat'])
    akv40 = u;      akv40 = squeeze(akv40(:,76));
    akv40v = v;     akv40v = squeeze(akv40v(:,76));

else
    rootdir = '\dataloc\';
    
    %Akv = 0.5e-3
    load([rootdir,'Akv_05\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km50g099994sx200it2000000.mat']);
    akv05 = u;      akv05 = squeeze(akv05(:,76));
    akv05v = v;     akv05v = squeeze(akv05v(:,76));
    
    %Akv = 0.75e-3
    load([rootdir,'Akv_075\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km75g099994sx200it2000000.mat']);
    akv075 = u;     akv075 = squeeze(akv075(:,76));
    akv075v = v;    akv075v = squeeze(akv075v(:,76));
    
    %Akv = 1.0e-3
    load([rootdir,'Akv_1\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km100g099994sx200it2000000.mat']);
    akv10 = u;      akv10 = squeeze(akv10(:,76));
    akv10v = v;     akv10v = squeeze(akv10v(:,76));
    
    %Akv = 1.25e-3
    load([rootdir,'Akv_125\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km125g10001sx200it2000000.mat']);
    akv125 = u;     akv125 = squeeze(akv125(:,76));
    akv125v = v;    akv125v = squeeze(akv125v(:,76));
    
    %Akv = 1.5e-3
    load([rootdir,'Akv_15\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km150g1sx200it2000000.mat'])
    akv15 = u;      akv15 = squeeze(akv15(:,76));
    akv15v = v;     akv15v = squeeze(akv15v(:,76));
    
    %Akv = 1.75e-3
    load([rootdir,'Akv_175\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km175g1sx200it2000000.mat']);
    akv17 = u;      akv17 = squeeze(akv17(:,76));
    akv17v = v;     akv17v = squeeze(akv17v(:,76));
    
    %Akv = 2.0e-3
    load([rootdir,'Akv_2\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km200g1sx200it2000000.mat'])
    akv20 = u;      akv20 = squeeze(akv20(:,76));
    akv20v = v;     akv20v = squeeze(akv20v(:,76));

    %Akv = 2.25e-3
    load([rootdir,'Akv_225\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km225g15625sx200it2000000.mat']);
    akv225 = u;     akv225 = squeeze(akv225(:,76));
    akv225v = v;    akv225v = squeeze(akv225v(:,76));
    
    %Akv = 2.5e-3
    load([rootdir,'Akv_25\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km250g17361sx200it2000000.mat'])
    akv25 = u;      akv25 = squeeze(akv25(:,76));
    akv25v = v;     akv25v = squeeze(akv25v(:,76));
    
    %Akv = 2.75e-3
    load([rootdir,'Akv_275\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km275g19097sx200it2000000.mat']);
    akv27 = u;      akv27 = squeeze(akv27(:,76));
    akv27v = v;     akv27v = squeeze(akv27v(:,76));
    
    %Akv = 3.0e-3
    load([rootdir,'Akv_3\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km300g20833sx200it2000000.mat'])
    akv30 = u;      akv30 = squeeze(akv30(:,76));
    akv30v = v;     akv30v = squeeze(akv30v(:,76));

    %Akv = 3.25e-3
    load([rootdir,'Akv_325\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km325g22569sx200it2000000.mat']);
    akv325 = u;     akv325 = squeeze(akv325(:,76));
    akv325v = v;    akv325v = squeeze(akv325v(:,76));
    
    %Akv = 3.5e-3
    load([rootdir,'Akv_35\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km350g24306sx200it2000000.mat'])
    akv35 = u;      akv35 = squeeze(akv35(:,76));
    akv35v = v;     akv35v = squeeze(akv35v(:,76));
    
    %Akv = 3.75e-3
    load([rootdir,'Akv_375\ModelNumerics_NonRot\rk_suvw_nx151B10000h12Km375g26042sx200it2000000.mat']);
    akv37 = u;      akv37 = squeeze(akv37(:,76));
    akv37v = v;     akv37v = squeeze(akv37v(:,76));
    
    %Akv = 4.0e-3
    load([rootdir,'Akv_4\ModelNumerics\rk_suvw_nx151B10000h12Km400g27778sx200it2000000.mat'])
    akv40 = u;      akv40 = squeeze(akv40(:,76));
    akv40v = v;     akv40v = squeeze(akv40v(:,76));
end

%% - Concantenate Date

blank = zeros(size(akv20,1),size(akv20,2)).*nan;

uarray = cat(2,akv05,akv075,akv10,akv125,akv15,akv17,akv20,akv225,akv25,akv27,akv30,akv325,akv35,akv37,akv40);
varray = cat(2,akv05v,akv075v,akv10v,akv125v,akv15v,akv17v,akv20v,akv225v,akv25v,akv27v,akv30v,akv325v,akv35v,akv37v,akv40v);

x1D = 0.0005:0.00025:0.004;
x2D = ones(size(uarray,1),size(uarray,2)).*x1D;
x2Dv = ones(size(varray,1),size(varray,2)).*x1D;

%% - Plot Line

if coriolis == 1
    vplot = zeros(1,length(x1D));
    test=0;
    for i = 1:length(x1D)
        ucut = squeeze(uarray(47:106,i));
        vcut = squeeze(varray(47:106,i));
        indi=find(diff(sign(vcut)));
        if length(indi) > 1
            indi = 47;
            test = test+1;
        end
        vplot(1,i) = squeeze(ucut(indi,1));
    end
else
    vplot = squeeze(uarray(76,:));
end

coefficients = polyfit(x1D,vplot,3);
xFit = linspace(min(x1D),max(x1D),1000);
yFit = polyval(coefficients,xFit);

figure
plot(x1D,vplot,'k-')
hold on
plot(x1D,vplot,'r.')
plot(xFit,yFit,'c-')
xlim([0.0005 0.004])
xlabel('Vertical Eddy Viscosity (m^2 s^{-1})')
%ylim([-0.5 0.5])
ylabel('Along-Channel Velocity (m s^{-1})')
title('Center Channel Along-Channel Velocity as a function of Akv')
grid on
if isv == 1
    print([rootdir,'CenterSurfacePlot'],'-dpng');
end

%% - Plot pcolor

widthplotr = ones(size(uarray,1),size(uarray,2)).*(((0:(ny-1)).*dy)');
acrossdist = ((ones(size(uarray,1),size(uarray,2)).*widthplotr)-5000)./1000;
widthplotv = ones(size(varray,1),size(varray,2)).*(((0:ny).*dy)');
acrossdisv = ((ones(size(varray,1),size(varray,2)).*widthplotv)-5000)./1000;

cmap = cmocean('balance');

figure
set(gcf,'Position',[150  150  800  600])
pcolor(x2D,acrossdist,-uarray)
shading('flat')
hold on
conval = -1:0.025:1;
[C,h] = contour(x2Dv,acrossdisv,varray,conval,'k-','ShowText','on');
clabel(C,h,'FontSize',15,'Color','black')
xlim([0.0005 0.004])
xlabel('Vertical Eddy Viscosity (m^2 s^{-1})','FontSize',16)
%ylim([-5 5])
ylabel('Across-Channel Distance (km)','FontSize',16)
ax = gca; 
ax.FontSize = 16;
title('Surface Along-Channel Velocity as a function of Akv')
grid on
c1 = colorbar;
colormap(cmap)
caxis([-0.06 0.06]);
c1.Label.String = 'Along-Channel Velocity (m s^{-1})';
if isv == 1
    print([rootdir,'FullSurfacePlot'],'-dpng');
end


