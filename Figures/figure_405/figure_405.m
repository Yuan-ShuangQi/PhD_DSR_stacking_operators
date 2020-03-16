%% Accuracy of the DSR-PS stacking operator

%% Introduction
% 
% *Author*: Abakumov Ivan
%
% *Publication date*: 1st April 2016


%% Define working folder, add links to Library and SeisLab 

clear; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;
current_folder = pwd;

%% Define model parameters (See table 4.1): 

alpha = pi/6;               % radian 
Rnip = 0.5;                 % km 
Rn   = 1.0;                 % km
Vp   = 2.5;                 % km/s 
Vs   = 1.8;                 % km/s   
modelPS = [alpha, Rnip, Rn, Vp, Vs]; 
modelPP = [alpha, Rnip, Rn, Vp, Vp]; 

gamma = Vp/Vs; 
sigma = (gamma-1)/(gamma+1);
upsilon = 2/(gamma+1);

%% Set offset and midpoint displacement

m = -0.5:0.010:0.5; 
h = -0.0:0.010:2.0;

[M,H]=meshgrid(m,h); 

%make gamma-CMP coordinates
tH = upsilon*H; 
tM = M + sigma*H; 

%% Part I: Calculate traveltimes of PP and PS waves

[TePP,~] = Get_traveltime_2D_exact(M,H,modelPP); 
[TePS,~] = Get_traveltime_2D_exact(M,H,modelPS); 
[T6PS,~] = Get_traveltime_2D_sixth_order(M, H, modelPS);

%% Part II: Find traveltime approximation for PP and PS waves

% PP approximations
T_CRS     = Get_traveltime_2D_CRS(M, H, modelPP);
T_MF      = Get_traveltime_2D_MF(M, H, modelPP);
T_iCRS    = Get_traveltime_2D_iCRS(M, H, modelPP);
T_nCRS    = Get_traveltime_2D_nCRS(M, H, modelPP);
T_DSR_PP  = Get_traveltime_2D_DSR_PS(M, H, modelPP);

% PS approximations
T_DSR_PS  = Get_traveltime_2D_DSR_PS(M, H, modelPS);
T_CRS_PS  = Get_traveltime_2D_CRS_PS(M, H, modelPS);
T_nCRS_PS = Get_traveltime_2D_nCRS_PS(M, H, modelPS);


%% Part III: Compare PP approximations

% Plot results
figure(1)
%
subplot(2,3,1)
data = (T_DSR_PP'-TePP')./TePP'*100;
imagesc(h,m',flipud(data))
title('DSR')
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
colorbar
colormap('Jet'); 
caxis([-5 5])
%
subplot(2,3,2)
data = (T_CRS'-TePP')./TePP'*100;
imagesc(h,m',flipud(data))
title('CRS')
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
colorbar
colormap('Jet'); 
caxis([-5 5])
%
subplot(2,3,4)
data = (T_MF'-TePP')./TePP'*100;
imagesc(h,m',flipud(data))
title('MF')
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
colorbar
colormap('Jet'); 
caxis([-5 5])
%
subplot(2,3,5)
data = (T_iCRS'-TePP')./TePP'*100;
imagesc(h,m',flipud(data))
title('i-CRS')
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
colorbar
colormap('Jet'); 
caxis([-5 5])
%
subplot(2,3,6)
data = (T_nCRS'-TePP')./TePP'*100;
imagesc(h,m',flipud(data))
title('n-CRS')
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
colorbar
colormap('Jet'); 
caxis([-5 5])

%% Part IV Compare PS approximations (old figure)

ind = (m==0); 
texact = TePP(:, ind); 
tMF    = T_MF(:, ind); 
tCRS   = T_CRS(:, ind); 
tnCRS  = T_nCRS(:, ind); 
tDSR   = T_DSR_PP(:, ind); 


figure(2)
subplot(2,2,1)
plot(h, texact, 'black'); 
hold on 
plot(h, tMF, 'g'); 
plot(h, tCRS, 'blue'); 
plot(h, tnCRS, 'm');
plot(h, tDSR, 'r');
legend('Exact', 'MF', 'CRS', 'n-CRS', 'DSR');
axis([0, 1.35, 0.35, 1.35]); 
xlabel('Offset [km]')
ylabel('Traveltime [sec]')
title('Exact traveltime, MF, CRS n-CRS and DSR approximations')

subplot(2,2,3)
plot(h, (tMF-texact)./texact*100, 'g'); 
hold on
plot(h, (tCRS-texact)./texact*100, 'blue'); 
plot(h, (tnCRS-texact)./texact*100, 'm');
plot(h, (tDSR-texact)./texact*100, 'r');
legend('MF', 'CRS', 'n-CRS', 'DSR');
axis([0, 1.35, -10, 10]); 
xlabel('Offset [km]')
ylabel('Error, %')
title('MF, CRS, n-CRS and DSR errors')

subplot(2,2,2)
data = (T_DSR_PS'-TePS')./TePS'*100;
imagesc(h,m,data)
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
title('Relative error of DRS-PS approximation, PS wave (%)')
colorbar
colormap('Jet'); 
caxis([-1 5])

subplot(2,2,4)
filter = (abs(tM)<0.2); 
data = (T_DSR_PS'-TePS')./TePS'*100.*filter'; 
data = data + 100*(filter'-1); 
imagesc(h,m,data)
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
title('Relative error of DRS-PS approximation, PS wave (%)')
colorbar
colormap('Jet'); 
caxis([-1 1])

%% Part V Compare PS approximations (new figure)

ind = (m==0); 
texact  = TePS(:, ind); 
tCRSPS  = T_CRS_PS(:, ind); 
tnCRSPS = T_nCRS_PS(:, ind); 
tDSRPS  = T_DSR_PS(:, ind); 

figure(3)
subplot(2,2,1)
plot(h, texact, '-black', 'LineWidth',2); 
hold on 
plot(h, tCRSPS, '-blue', 'LineWidth',2); 
plot(h, tnCRSPS, '-g', 'LineWidth',2);
plot(h, tDSRPS, '-r', 'LineWidth',2);
legend('Exact', 'CRS-PS', 'n-CRS-PS', 'DSR-PS','Location','NorthWest');
axis([0, 1.35, 0.35, 1.35]); 
xlabel('Offset [km]')
ylabel('Traveltime [sec]')
title('Exact traveltime, CRS-PS, n-CRS-PS and DSR-PS approximations')

subplot(2,2,3)
plot(h, (tCRSPS-texact)./texact*100, '-blue', 'LineWidth',2); 
hold on
plot(h, (tnCRSPS-texact)./texact*100, '-g', 'LineWidth',2);
plot(h, (tDSRPS-texact)./texact*100, '-r', 'LineWidth',2);
plot(h, (texact-texact)./texact*100, '--black', 'LineWidth',1); 
legend('CRS-PS', 'n-CRS-PS', 'DSR-PS','Location','NorthWest');
axis([0, 1.35, -10, 10]); 
xlabel('Offset [km]')
ylabel('Error, %')
title('CRS-PS, n-CRS-PS and DSR-PS errors')

subplot(2,2,2)
data = (T_DSR_PS'-TePS')./TePS'*100;
imagesc(h,m',flipud(data))
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
title('Relative error of DSR-PS approximation (%)')
colorbar
colormap('Jet'); 
caxis([-5 5])
hold on 
contour(h,m',flipud(data),'color','k','ShowText','on')

subplot(2,2,4)
data = (T_nCRS_PS'-TePS')./TePS'*100;
imagesc(h,m',flipud(data))
xlabel('Offset [km]'); 
ylabel('Mid-point [km]');
title('Relative error of n-CRS-PS approximation (%)')
colorbar
colormap('Jet'); 
caxis([-5 5])
hold on 
contour(h,m',flipud(data),'color','k','ShowText','on')