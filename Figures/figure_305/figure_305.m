%% Figure 305
% Check 3D i-CRS TIA/LIA algorithms
% 
% *Author*: Abakumov Ivan
% *Publication date*: 31st August 2016

%% Define working folder, add links to Library and SeisLab 

clear; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;

%% Introduction

% Note:
% Accuracy of traveltime 10e-12
% Accuracy of attributes 10e-10

for model = 0:0
    for acquisition = 3:3

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

     [ t0, w3, M3, N3 ] = Get_model_stacking_parameters( X0, model );
      
     CRS_param.x0 = X0; 
     CRS_param.t0 = t0; 
     CRS_param.v0 = v0; 
     CRS_param.w  = w3; 
     CRS_param.M  = M3; 
     CRS_param.N  = N3; 
      
     save([mlibfolder '/CRS/models/model_' num2str(model) '_CRS_param.mat'], 'CRS_param'); 
 
     % either take (x,y) components and find wavefield attributes:
     % (reality)
     w = w3(1:2, 1); 
     N = N3(1:2, 1:2);
     M = M3(1:2, 1:2);  
     %[ alpha, beta, KNIP, KN ] = my_A2P(v0, w, M, N);
         
     % or find wavefield attributes and make from them stacking parameters:
     % (theory)
     % [ alpha, beta, KNIP3, KN3 ] = my_A3P(v0, w3, M3, N3);
     %  KNIP = KNIP3(1:2,1:2); 
     %  KN   = KN3  (1:2,1:2);
     %  [ w, M, N ] = my_P2A(v0, alpha, beta, KNIP, KN);
    
    %% Exact traveltimes
    
    % download exact traveltimes
    % tti_ex = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
    
    % or calculation exact traveltimes again

     tic
     tti_ex = Get_model_exact_traveltime(Xs, Xg, model);
     toc
     save([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat'], 'tti_ex'); 
       
    %% i-CRS approximations
    
     HH = (Xg(1:2, :) - Xs(1:2,:))/2;
     MM = (Xg(1:2, :) + Xs(1:2,:))/2;
     MM(1,:) = MM(1,:) - X0(1); 
     MM(2,:) = MM(2,:) - X0(2);
     
     %
     tic
     tti_icrs_par_LIA = Get_traveltime_3D_iCRS_par_LIA(MM, HH, t0, w, M, N, 20);
     ctime.icrs_par_LIA = toc; 
     %
     tic
     tti_icrs_el_LIA = Get_traveltime_3D_iCRS_el_LIA(MM, HH, t0, v0, w, M, N, 20);
     ctime.icrs_el_LIA = toc; 
     %
     tic
     tti_icrs_el_TIA = Get_traveltime_3D_iCRS_el_TIA(MM, HH, t0, v0, w, M, N, 20);
     ctime.icrs_el_LTA = toc; 

    %% RMS errors

    for iter = 1:20
        err_rms_icrs_par_LIA(iter) = sqrt(sum(((tti_icrs_par_LIA(iter,:) - tti_ex)./tti_ex).^2)/100)*100;
        err_rms_icrs_el_LIA(iter)  = sqrt(sum(((tti_icrs_el_LIA(iter,:)  - tti_ex)./tti_ex).^2)/100)*100;
        err_rms_icrs_el_TIA(iter)  = sqrt(sum(((tti_icrs_el_TIA(iter,:)  - tti_ex)./tti_ex).^2)/100)*100;
    end

    %% Plot figure 305
    subplot(1,2,1)
    plot(Xs(1,:), Xs(2,:), 'rV')
    hold on 
    plot(Xg(1,:), Xg(2,:), 'g^')
    plot(X0(1), X0(2), '*black')
    legend('Shot', 'Receiver', 'Central point')
    for a = 0:pi/6:pi
        hold on
        x1 = (-1000:10:1000)*cos(a) + X0(1);
        x2 = (-1000:10:1000)*sin(a) + X0(2);
        plot(x1,x2,'Color',[0.75 0.75 0.75])
    end
    plot(1000*sin(0:pi/100:2*pi)+X0(1), 1000*cos(0:pi/100:2*pi)+X0(2),'--black') 
    plot(X0(1), X0(2), '*black')
    text(3100,2250,'1000')
    xlabel('Inline  [m]'); 
    ylabel('Crossline [m]'); 
    axis('equal')
    axis([500, 3500, 500, 3500])

    subplot(1,2,2); 
    
    semilogy(err_rms_icrs_el_LIA, '--black')
    hold on 
    semilogy(err_rms_icrs_el_TIA, '--red')
    legend('LIA', 'TIA')
    xlabel('Iteration number'); 
    ylabel('Error %'); 
    axis('square')
 
    end
end

