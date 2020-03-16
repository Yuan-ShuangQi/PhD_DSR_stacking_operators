%% Figure 314
% Calculate computational time of CRS / DSR / nCRS and i-CRS stacking operators
% 
% *Author*: Abakumov Ivan
% *Publication date*: 17st October 2016

%% Define working folder, add links to Library and SeisLab 

clear; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;

%% Introduction

% Note:
% Accuracy of traveltime 10e-12
% Accuracy of attributes 10e-10

for model = 63:63
    acquisition = 2;

    %% Get model parameters
    
    Get_model_parameters; 
    
    %% Get acquisition geometry

    Get_model_acquisition_geometry; 

    %% Plot acquisition geometry
    figure(1)
    [XX, YY] = meshgrid(G.xx, G.yy);
    [ZZ, ind] = Get_model_surface(XX,YY,model);
    surf(XX,YY,ZZ); 
    hold on
    plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
    plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
    plot3(X0(1),X0(2),X0(3), 'b*');
    axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
    xlabel('Inline  [m]'); 
    ylabel('Crossline [m]'); 
    zlabel('Depth  [m]'); 
    view(3); 
    set(gca, 'ZDir', 'reverse')
    
    %% Find stacking parameters

     CRS_param = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_CRS_param.mat']); 

     X0 = CRS_param.x0; 
     t0 = CRS_param.t0; 
     v0 = CRS_param.v0; 
     w3 = CRS_param.w; 
     M3 = CRS_param.M; 
     N3 = CRS_param.N; 
      
     % either take (x,y) components and find wavefield attributes:
     % (reality)
     w = w3(1:2, 1); 
     N = N3(1:2, 1:2);
     M = M3(1:2, 1:2);  
     %[ alpha, beta, KNIP, KN ] = my_A2P(v0, w, M, N);
         
     % or find wavefield attributes and make from them stacking parameters:
     % (theory)
%      [ alpha, beta, KNIP3, KN3 ] = my_A3P(v0, w3, M3, N3);
%       KNIP = KNIP3(1:2,1:2); 
%       KN   = KN3  (1:2,1:2);
%       [ w, M, N ] = my_P2A(v0, alpha, beta, KNIP, KN);
%     
    %% Exact traveltimes
    
    % download exact traveltimes
    tti_ex = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
    
    % or calculation exact traveltimes again

  %    tic
  %    tti_ex = Get_model_exact_traveltime(Xs, Xg, model);
  %    toc
  %    save([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat'], 'tti_ex'); 
       
    %% Traveltime approximations
    
    HH = (Xg(1:2, :) - Xs(1:2,:))/2;
    MM = (Xg(1:2, :) + Xs(1:2,:))/2;
    MM(1,:) = MM(1,:) - X0(1); 
    MM(2,:) = MM(2,:) - X0(2);
     
    %% temp
    
    ctime.crs = 0;
    ctime.dsr = 0;
    ctime.icrs_el_LIA = 0;
    ctime.icrs_el_TIA = 0;
    ctime.icrs_par_LIA = 0;
    ctime.ncrs = 0;
     
     %%
    tic
    tti_crs = Get_traveltime_3D_CRS (MM, HH, t0, w, M, N);
    ctime.crs = toc; 
    %
    tic
    tti_dsr = Get_traveltime_3D_DSR (MM, HH, t0, w, M, N);
    ctime.dsr = toc; 
    %
    tic
    tti_ncrs = Get_traveltime_3D_nCRS(MM, HH, t0, w, M, N);
    ctime.ncrs = toc; 
    %
    tic
    tti_icrs_par_LIA = Get_traveltime_3D_iCRS_par_LIA(MM, HH, t0, w, M, N, 1);
    ctime.icrs_par_LIA = toc; 
    %
    tic
    tti_icrs_el_LIA = Get_traveltime_3D_iCRS_el_LIA(MM, HH, t0, v0, w, M, N, 1);
    ctime.icrs_el_LIA = toc; 
    %
    tic
    tti_icrs_el_TIA = Get_traveltime_3D_iCRS_el_TIA(MM, HH, t0, v0, w, M, N, 1);
    ctime.icrs_el_TIA = toc;    
end

