%% Figure 301
% Discover accuracy of the simplified model
% 
% *Author*: Abakumov Ivan
% *Publication date*: 31st August 2016

%% Define working folder, add links to Library and SeisLab 

clear; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;


%% Plot results for model 63

 
% Download exact traveltimes
model = 63; 
acquisition = 4; 
tti_ex_CMP = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
acquisition = 5; 
tti_ex_ZO = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
    
% Download exact traveltimes in the simplified model with paraboloid
model = 163;
acquisition = 4; 
tti_SM_par_CMP = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
acquisition = 5; 
tti_SM_par_ZO = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
     
% Download exact traveltimes in the simplified model with ellipsoid
model = 263;
acquisition = 4; 
tti_SM_el_CMP = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
acquisition = 5;
tti_SM_el_ZO = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
 
% Calculate CRS approximations
model = 63;
CRS_param = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_CRS_param.mat']); 
t0 = CRS_param.t0; 
w3  = CRS_param.w; 
M3  = CRS_param.M; 
N3  = CRS_param.N;
w   = w3(1:2); 
M   = M3(1:2,1:2); 
N   = N3(1:2,1:2); 
 
acquisition = 4; 
Get_model_acquisition_geometry; 
HH_CMP = (Xg(1:2, :) - Xs(1:2,:))/2;
MM_CMP = (Xg(1:2, :) + Xs(1:2,:))/2;
MM_CMP(1,:) = MM_CMP(1,:) - X0(1); 
MM_CMP(2,:) = MM_CMP(2,:) - X0(2);
of_CMP = offset; 
az_CMP = azimuth; 
tti_crs_CMP  = Get_traveltime_3D_CRS (MM_CMP, HH_CMP, t0, w, M, N);
     
acquisition = 5; 
Get_model_acquisition_geometry; 
HH_ZO = (Xg(1:2, :) - Xs(1:2,:))/2;
MM_ZO = (Xg(1:2, :) + Xs(1:2,:))/2;
MM_ZO(1,:) = MM_ZO(1,:) - X0(1); 
MM_ZO(2,:) = MM_ZO(2,:) - X0(2);
of_ZO = offset; 
az_ZO = azimuth; 
tti_crs_ZO  = Get_traveltime_3D_CRS (MM_ZO, HH_ZO, t0, w, M, N);  
  
%%

dt_crs_CMP    = reshape((tti_crs_CMP   -tti_ex_CMP)./tti_ex_CMP*100,length(az_CMP),length(of_CMP));
dt_SM_el_CMP  = reshape((tti_SM_el_CMP -tti_ex_CMP)./tti_ex_CMP*100,length(az_CMP),length(of_CMP));
dt_SM_par_CMP = reshape((tti_SM_par_CMP-tti_ex_CMP)./tti_ex_CMP*100,length(az_CMP),length(of_CMP));

dt_crs_ZO    = reshape((tti_crs_ZO   -tti_ex_ZO)./tti_ex_ZO*100,length(az_ZO),length(of_ZO));
dt_SM_el_ZO  = reshape((tti_SM_el_ZO -tti_ex_ZO)./tti_ex_ZO*100,length(az_ZO),length(of_ZO));
dt_SM_par_ZO = reshape((tti_SM_par_ZO-tti_ex_ZO)./tti_ex_ZO*100,length(az_ZO),length(of_ZO));
 
[OF, AZ] = meshgrid(of_CMP,az_CMP);
XX = OF.*cos(AZ); 
YY = OF.*sin(AZ); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,1)
% h = polar(azimuth,2000*ones(size(azimuth)), '-black');
% hold on
% contourf(XX,YY,dt_crs_CMP,11);
% %
% for a = 0:pi/6:pi
%     hold on
%     x1 = (-2000:10:2000)*cos(a);
%     x2 = (-2000:10:2000)*sin(a);
%     plot(x1,x2,'Color',[0.75 0.75 0.75])
% end
% polar(azimuth,1000*ones(size(azimuth)),'--black');
% text(250,1100,'1000')
% text(-2000,2000,'a)')
% %
% c = colorbar('Ticks', [-1,-0.5 0, 0.5, 1]);
% c.Label.String = 'Error %';
% caxis([-1 1])
% xlabel('h_x')
% ylabel('h_y')
% title('CRS approximation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
h = polar(azimuth,2000*ones(size(azimuth)), '-black');
hold on
contourf(XX,YY,dt_SM_el_CMP,11);
for a = 0:pi/6:pi
    hold on
    x1 = (-2000:10:2000)*cos(a);
    x2 = (-2000:10:2000)*sin(a);
    plot(x1,x2,'Color',[0.75 0.75 0.75])
end
polar(azimuth,1000*ones(size(azimuth)),'--black');
text(250,1100,'1000')
text(-2000,2000,'a)')
c = colorbar('Ticks', [-1,-0.5 0, 0.5, 1]);
c.Label.String = 'Error %';
caxis([-1 1])
%colormap(makeColorMap([0 0 1],[1 1 1],[1 0 0],100));
xlabel('h_x [m]')
ylabel('h_y [m]')
title('CMP experiment, ellipsoidal reflector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
h = polar(azimuth,2000*ones(size(azimuth)), '-black');
hold on
contourf(XX,YY,dt_SM_par_CMP,11);
%
for a = 0:pi/6:pi
    hold on
    x1 = (-2000:10:2000)*cos(a);
    x2 = (-2000:10:2000)*sin(a);
    plot(x1,x2,'Color',[0.75 0.75 0.75])
end
polar(azimuth,1000*ones(size(azimuth)),'--black');
text(250,1100,'1000')
text(-2000,2000,'b)')
c = colorbar('Ticks', [-1,-0.5 0, 0.5, 1]);
c.Label.String = 'Error %';
caxis([-1 1])
xlabel('h_x [m]')
ylabel('h_y [m]')
title('CMP experiment, parabolic reflector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[OF, AZ] = meshgrid(of_ZO,az_ZO);
XX = OF.*cos(AZ); 
YY = OF.*sin(AZ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,2,4)
% h = polar(azimuth,1000*ones(size(azimuth)), '-black');
% hold on
% contourf(XX,YY,dt_crs_ZO,11);
% %
% for a = 0:pi/6:pi
%     hold on
%     x1 = (-1000:10:1000)*cos(a);
%     x2 = (-1000:10:1000)*sin(a);
%     plot(x1,x2,'Color',[0.75 0.75 0.75])
% end
% polar(azimuth,500*ones(size(azimuth)),'--black');
% text(125,550,'500')
% text(-1000,1000,'d)')
% c = colorbar('Ticks', [-2,-1 0, 1, 2]);
% c.Label.String = 'Error %';
% caxis([-2 2])
% xlabel('m_x')
% ylabel('m_y')
% title('CRS approximation')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
h = polar(azimuth,1000*ones(size(azimuth)), '-black');
hold on
contourf(XX,YY,dt_SM_el_ZO,11);
%
for a = 0:pi/6:pi
    hold on
    x1 = (-1000:10:1000)*cos(a);
    x2 = (-1000:10:1000)*sin(a);
    plot(x1,x2,'Color',[0.75 0.75 0.75])
end
polar(azimuth,500*ones(size(azimuth)),'--black');
text(125,550,'500')
text(-1000,1000,'c)')
c = colorbar('Ticks', [-2,-1 0, 1, 2]);
c.Label.String = 'Error %';
caxis([-2 2])
xlabel('m_x [m]')
ylabel('m_y [m]')
title('ZO experiment, ellipsoidal reflector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)
h = polar(azimuth,1000*ones(size(azimuth)), '-black');
hold on
contourf(XX,YY,dt_SM_par_ZO,11);
%
for a = 0:pi/6:pi
    hold on
    x1 = (-1000:10:1000)*cos(a);
    x2 = (-1000:10:1000)*sin(a);
    plot(x1,x2,'Color',[0.75 0.75 0.75])
end
polar(azimuth,500*ones(size(azimuth)),'--black');
text(125,550,'500')
text(-1000,1000,'d)')
colorbar;  
c = colorbar('Ticks', [-2,-1 0, 1, 2]);
c.Label.String = 'Error %';
caxis([-2 2])
xlabel('m_x [m]')
ylabel('m_y [m]')
title('ZO experiment, parabolic reflector')
    