%% figure p01
%
% *Author*: Abakumov Ivan
%
% *Publication date*: 18th August 2016

%% Define working folder, add links to Library and SeisLab 

clear all; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;

%% Read data
data = read_su_file('/home/zmaw/u250128/Desktop/CRS/SEG_C3_WA/3D_SEG_cdp_181903_182236.su');

header = get_segy_headers(data);

 
scalco = 1000;                            
scalel = 10000; 
  
header.sx = header.sx/scalco; 
header.sy = header.sy/scalco; 
header.rx = header.rx/scalco; 
header.ry = header.ry/scalco; 
header.sz = header.sz/scalel; 
header.rz = header.rz/scalel; 

change_scale = (header.sx < 3700); 
header.sx(change_scale) = header.sx(change_scale)*10; 
header.sy(change_scale) = header.sy(change_scale)*10; 
header.rx(change_scale) = header.rx(change_scale)*10; 
header.ry(change_scale) = header.ry(change_scale)*10; 

header.cdpx = (header.sx+header.rx)/2;
header.cdpy = (header.sy+header.ry)/2; 
header.hx = (header.rx - header.sx)/2; 
header.hy = (header.ry - header.sy)/2;

%%
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
DIS = (G.xx(103:436)-G.xx(103))/1000;
DEP = G.zz/1000;
pcolor(DIS,DEP, squeeze(velmod(400,103:436,:))'/1000)
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Depth, km','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6 8 10 12], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
c = colorbar('eastoutside', 'Ticks', [1,2,3,4]);
shading interp
set(gca, 'YDir', 'reverse')
c.Label.String = 'Velocity, km/s';
c.Label.Color = [51/255 51/255 153/255];
c.Color = [51/255 51/255 153/255];
axis('equal')
axis([min(DIS) max(DIS) 0 3.5]);
caxis([1 3.5]); 
colormap('jet')

%%
ind = header.cdpy == 4500; 

hold on
plot(header.ry(ind)/1000, 0, '*g')
plot(header.sy(ind)/1000, 0, '*r')

%%
seismic = data; 
seismic.traces  = data.traces(:,ind); 
seismic.headers = data.headers(:,ind); 
seismic.headers(7,:) = round(seismic.headers(7,:)/100)/10;
f = s_wplot(seismic,{'annotation','offset'},{'figure','new'},{'quality','draft'},{'title',''});
xlabel('Offset, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'xtick',[], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
title('');
hold on
tt = sqrt(1 + offset.^2/1500^2);
plot(1:17, 1000*tt, '--r', 'Linewidth',3.0)
tt = sqrt(1 + offset.^2/2000^2);
plot(1:17, 1000*tt, '--r', 'Linewidth',3.0)
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off';

