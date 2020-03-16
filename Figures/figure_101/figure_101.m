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

%% Upload original prestack data

sufilename = '/scratch/local1/ivan/CRS_data/PS_Muehldorf_2D_data/Original_Muehldorf_data/Muehldorf_1/MUEHL8542P_CDP.SEGY';
data = read_segy_file(sufilename);
header = get_segy_headers(data);
header.ry = zeros(size(header.sx));
     
% cdp every 15 meters
header.sx = header.sx - 13940; 
header.rx = header.rx - 13940; 
header.cdpx = (header.sx+header.rx)/2;


%% Plot all acquisition geometry
figure(21)
subplot(1,6,[1 2 3])
plot(header.sx, header.rx, '.', 'color', [0.5 0.5 0.5])
hold on
indCS = (header.sx == 775); 
plot(header.sx(indCS), header.rx(indCS), '.r')
indCR = (header.rx == 1000); 
plot(header.sx(indCR), header.rx(indCR), '.g')
indCMP = (header.cdp == 2927); 
plot(header.sx(indCMP), header.rx(indCMP), '.b')
axis([-1500 3000 -1500 3000]);
xlabel('Shot [m]');
ylabel('Receiver [m]')
title('Stacking chart'); 
legend('Multi-coverage seismic data', 'CS gather', 'CR gather','CMP gather','NorthWest');

%axis('equal')     
subplot(1,6,[4 5])
alltrack = 1:8766;
trackCMP = alltrack(indCMP);
tr = 0:0.002:3.998; 
for i = 1:length(trackCMP)
    plot(header.offset(trackCMP(i))+15*data.traces(:,trackCMP(i)), tr, 'black');
    hold on
end
h = header.offset(trackCMP);
t0 = 1.4; 
vNMO = 5500; 
t = sqrt(t0^2 + 4*h.^2/vNMO^2); 
hold on 
plot(h, t, '-red', 'LineWidth', 2); 
set(gca, 'YDir', 'reverse'); 
xlabel('Half-offset [m]');
ylabel('Traveltime [sec]')
title('CMP gather'); 



  
%% Upload stacked data

subplot(1,6,6)
Ssufilename = '/scratch/local1/ivan/CRS_data/PS_Muehldorf_2D_data/Original_Muehldorf_data/Muehldorf_1/MUEHL8542P_STA.SEGY';
Sdata = read_segy_file(Ssufilename);
Sheader = get_segy_headers(Sdata);
indSTA = (Sheader.cdp == 2927); 
plot(Sdata.traces(:,indSTA)/10, tr, 'black');
axis([-1 1 0 4])
set(gca, 'YDir', 'reverse'); 
xlabel('Amplitude');
ylabel('Traveltime [sec]')
title('Stacked trace'); 
