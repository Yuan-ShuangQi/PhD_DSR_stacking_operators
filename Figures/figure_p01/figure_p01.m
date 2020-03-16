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

%% 
cdpmin = 181903;
cdpmax = 182236;
tt = (0:624)*0.008; 
cc = (0:333)*0.02; 

%% Plot CMP stack
fig = gcf;
fig.Color = [1.0 1.0 243/255];
%fig.InvertHardcopy = 'off';



data_cmp = read_su_file('/scratch/local1/ivan/CRS_data/SEG_C3_WA_3D_CRS/product_SEG_C3WA_ffid_proc_red.su_323/cmp/cmpsa.stack.su'); 
headers_cmp = get_segy_headers(data_cmp); 
[cdp, ind] = sort(headers_cmp.cdp);
traces_cmp = data_cmp.traces(:,ind); 
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind = (ind1+ind2) == 2; 

subplot(1,3,1)
pcolor(cc, tt, traces_cmp(:,ind))
caxis([-1 1])
colormap('gray')
title('CMP stack','FontSize',25, 'Color',[51/255 51/255 153/255])
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
axis([min(cc) max(cc) 0.16 4.5])
set(gca, 'YDir', 'reverse')
shading interp

%% Plot CRS stack

data_crs = read_su_file('/scratch/local1/ivan/CRS_data/SEG_C3_WA_3D_CRS/product_SEG_C3WA_ffid_proc_red.su_323/crs/crssa.stack.su');
headers_crs = get_segy_headers(data_crs); 
[cdp, ind] = sort(headers_crs.cdp);
traces_crs = data_crs.traces(:,ind); 
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind = (ind1+ind2) == 2; 

subplot(1,3,2)
pcolor(cc, tt, traces_crs(:,ind))
caxis([-0.4 0.4])
colormap('gray')
title('CRS stack','FontSize',25, 'Color',[51/255 51/255 153/255])
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
axis([min(cc) max(cc) 0.16 4.5])
set(gca, 'YDir', 'reverse')
shading interp


%% Plot i-CRS stack

data_icrs = read_segy_file('/home/zmaw/u250128/Desktop/CRS/SEG_C3_WA/product/CRS/3D_SEG_cdp_181903_182236_crs.stack.sgy');

subplot(1,3,3)
pcolor(cc, tt, data_icrs.traces)
caxis([-1 1])
colormap('gray')
title('i-CRS stack','FontSize',25, 'Color',[51/255 51/255 153/255])
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
axis([min(cc) max(cc) 0.16 4.5])
set(gca, 'YDir', 'reverse')
shading interp