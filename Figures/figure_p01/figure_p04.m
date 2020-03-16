%% Figure 101
%
% *Author*: Abakumov Ivan
%
% *Publication date*: 14th April 2015

%% Introduction
% Figure 101 with Muehldorf data

%% Define working folder, add links to Library and SeisLab 
% Variable _ivnfolder_ defines path to the project's main folder. If this folder
% is moved, path _ivnfolder_ must be updated. _MLIB_ is a library with MATLAB
% and C++ scripts. Script |addmypath.m| adds specified folders to MATLAB
% search path. Script |definefolders.m| defines pathes to folders inside main folder. 

clear all; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;
current_folder = pwd;

%%
cdpmin = 181903;
cdpmax = 182236;
cdpcmp = 182150;
tt = (0:624)*0.008; 
cc = (0:333)*0.02; 


%% Upload CMP stacked data

data_cmp = read_su_file('/scratch/local1/ivan/CRS_data/SEG_C3_WA_3D_CRS/product_SEG_C3WA_ffid_proc_red.su_323/cmp/cmpsa.stack.su'); 
headers_cmp = get_segy_headers(data_cmp); 
[cdp, ind] = sort(headers_cmp.cdp);
traces_cmp = data_cmp.traces(:,ind); 


%% Plot CMP gather
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

subplot(1,5,[1 2])
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind = (ind1+ind2) == 2; 

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
colormap('gray')
title('CRS stack','FontSize',25, 'Color',[51/255 51/255 153/255])
xlabel('Distance, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 2 4 6], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
axis([min(cc) max(cc) 0.16 4.5])
set(gca, 'YDir', 'reverse')
shading interp







