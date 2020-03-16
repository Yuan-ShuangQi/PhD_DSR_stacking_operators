%% Figure 101
%
% *Author*: Abakumov Ivan
%
% *Publication date*: 14th January 2017

%% Introduction
% Figure 101 with Muehldorf data

%% Add MLIB folder

clear all; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;

%%
cdpmin = 181903;
cdpmax = 182236;
cdpcmp = 182150;
tt = (0:624)*0.008; 
cc = (0:333)*0.02; 

segyvelfile = '/scratch/local1/ivan/CRS_data/SEG_C3_WA_data/SEG_C3NA_Velocity.sgy';
veldata = read_segy_file(segyvelfile); 

inline = veldata.headers(6, :)/200+1;    % (20 m x 10 scaling factor)
xline  = veldata.headers(7, :)/200+1; 
    
% Make G-file
G=GridClass;
G.x0 = 0;    G.nx = 676;      G.dx = 20;   
G.y0 = 0;    G.ny = 676;      G.dy = 20; 
G.z0 = 0;    G.nz = 201;      G.dz = 20; 
G.t0 = 0;    G.nt = 625;      G.dt = 0.008;

G.gridInfo;
G.setGrid;
Gold = oldGrid(G);
    
% Make velocity cube
velmod = zeros(G.nx,G.ny,G.nz); 
for i=1:length(inline)
    velmod(inline(i), xline(i), :) = veldata.traces(:, i); 
end

%% Plot reference In-Line 401 (X==8000m)

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

subplot(1,5,[1 2 3])
DIS = (G.xx(103:436)-G.xx(103))/1000;
DEP = G.zz/1000;
pcolor(DIS,DEP, squeeze(velmod(400,103:436,:))'/1000)
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Depth, km','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6 8 10 12], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
c = colorbar('southoutside', 'Ticks', [1,2,3,4]);
shading interp
set(gca, 'YDir', 'reverse')
c.Label.String = 'Velocity, km/s';
c.Label.Color = [51/255 51/255 153/255];
c.Color = [51/255 51/255 153/255];
axis('equal')
axis([min(DIS) max(DIS) 0 3.5]);
caxis([1 3.5]); 
colormap('jet')








%% Upload CRS stack 

data_crs = read_su_file('/scratch/local1/ivan/CRS_data/SEG_C3_WA_3D_CRS/product_SEG_C3WA_ffid_proc_red.su_323/crs/crssa.stack.su');
headers_crs = get_segy_headers(data_crs); 
[cdp, ind] = sort(headers_crs.cdp);
traces_crs = data_crs.traces(:,ind); 
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind = (ind1+ind2) == 2; 

subplot(1,5,[4 5])
pcolor(cc, tt, traces_crs(:,ind))
caxis([-0.4 0.4])
%colormap('gray')
title('CRS stack','FontSize',25, 'Color',[51/255 51/255 153/255])
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
axis([min(cc) max(cc) 0.16 4.5])
set(gca, 'YDir', 'reverse')
shading interp







