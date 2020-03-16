%% Figure 307 
% Plot models 
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

acquisition = 3;


subplot(2,3,1)
model = 11;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX,YY,ZZ,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
plot3(X0(1),X0(2),X0(3), 'b*');
axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
xlabel('Inline  [m]'); 
ylabel('Crossline [m]'); 
zlabel('Depth  [m]'); 
title('Flat reflector')
view(3); 
set(gca, 'ZDir', 'reverse')
L = light;



subplot(2,3,2)
model = 21;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX,YY,ZZ,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
plot3(X0(1),X0(2),X0(3), 'b*');
axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
xlabel('Inline  [m]'); 
ylabel('Crossline [m]'); 
zlabel('Depth  [m]'); 
title('Plane dipping reflector')
view(3); 
set(gca, 'ZDir', 'reverse')
L = light;

 
subplot(2,3,3)
model = 31;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX,YY,ZZ] = ellipsoid(xc,yc,zc,80,80,80,50) ;
surf(XX,YY,ZZ,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
plot3(X0(1),X0(2),X0(3), 'b*');
axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
xlabel('Inline  [m]'); 
ylabel('Crossline [m]'); 
zlabel('Depth  [m]'); 
title('Point diffractor (R = 10 m)'); 
view(3); 
set(gca, 'ZDir', 'reverse')
L = light;

subplot(2,3,4)
model = 41;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX,YY,ZZ,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
plot3(X0(1),X0(2),X0(3), 'b*');
axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
xlabel('Inline  [m]'); 
ylabel('Crossline [m]'); 
zlabel('Depth  [m]'); 
title('Sphere (R = 1 km)'); 
view(3); 
set(gca, 'ZDir', 'reverse')
L = light;

subplot(2,3,5)
model = 51;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX,YY,ZZ,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
plot3(X0(1),X0(2),X0(3), 'b*');
axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
xlabel('Inline  [m]'); 
ylabel('Crossline [m]'); 
zlabel('Depth  [m]'); 
title('Ellipsoid'); 
view(3); 
set(gca, 'ZDir', 'reverse')
L = light;
 
subplot(2,3,6)
model = 61;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX,YY,ZZ,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:),Xs(2,:),Xs(3,:), 'rv');
plot3(Xg(1,:),Xg(2,:),Xg(3,:), 'g^');
plot3(X0(1),X0(2),X0(3), 'b*');
axis([G.x0, G.mx, G.y0, G.my, G.z0, G.mz]);
xlabel('Inline  [m]'); 
ylabel('Crossline [m]'); 
zlabel('Depth  [m]'); 
title('Complex surface'); 
view(3); 
set(gca, 'ZDir', 'reverse')
L = light;



%% Plot figure for presentation

acquisition = 3;
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

%-----------------------------------
subplot(2,3,1)
model = 11;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX/1000,YY/1000,ZZ/1000,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:)/1000,Xs(2,:)/1000,Xs(3,:)/1000, 'rv', 'Markersize', 10);
plot3(Xg(1,:)/1000,Xg(2,:)/1000,Xg(3,:)/1000, 'g^', 'Markersize', 10);
plot3(X0(1)/1000,X0(2)/1000,X0(3)/1000, 'b*', 'Markersize', 10);
axis([0 4 0 4 0 2])
xlabel('X  [km]', 'Fontsize', 16); 
ylabel('Y [km]','Fontsize', 16); 
zlabel('Depth  [km]','Fontsize', 16); 
title('Flat reflector','Fontsize', 28)
set(gca,'xtick',[0 2 ])
set(gca,'ytick',[0 2 4])
set(gca,'ztick',[0 1 2])
view(3); 
set(gca, 'ZDir', 'reverse')
set(gca, 'FontSize', 16); 
camlight('headlight') 
lighting gouraud

%-----------------------------------
subplot(2,3,2)
model = 21;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX/1000,YY/1000,ZZ/1000,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:)/1000,Xs(2,:)/1000,Xs(3,:)/1000, 'rv', 'Markersize', 10);
plot3(Xg(1,:)/1000,Xg(2,:)/1000,Xg(3,:)/1000, 'g^', 'Markersize', 10);
plot3(X0(1)/1000,X0(2)/1000,X0(3)/1000, 'b*', 'Markersize', 10);
axis([0 4 0 4 0 2])
xlabel('X  [km]', 'Fontsize', 16); 
ylabel('Y [km]','Fontsize', 16); 
zlabel('Depth  [km]','Fontsize', 16); 
title('Plane dipping reflector', 'Fontsize', 28)
set(gca,'xtick',[0 2 ])
set(gca,'ytick',[0 2 4])
set(gca,'ztick',[0 1 2])
view(3); 
set(gca, 'ZDir', 'reverse')
set(gca, 'FontSize', 16); 
camlight('headlight') 
lighting gouraud

%-----------------------------------
subplot(2,3,3)
model = 31;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX,YY,ZZ] = ellipsoid(xc,yc,zc,80,80,80,50);
surf(XX/1000,YY/1000,ZZ/1000,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:)/1000,Xs(2,:)/1000,Xs(3,:)/1000, 'rv', 'Markersize', 10);
plot3(Xg(1,:)/1000,Xg(2,:)/1000,Xg(3,:)/1000, 'g^', 'Markersize', 10);
plot3(X0(1)/1000,X0(2)/1000,X0(3)/1000, 'b*', 'Markersize', 10);
axis([0 4 0 4 0 2])
xlabel('X  [km]', 'Fontsize', 16); 
ylabel('Y [km]','Fontsize', 16); 
zlabel('Depth  [km]','Fontsize', 16); 
title('Point diffractor (R = 10 m)', 'Fontsize', 28); 
set(gca,'xtick',[0 2 ])
set(gca,'ytick',[0 2 4])
set(gca,'ztick',[0 1 2])
view(3); 
set(gca, 'ZDir', 'reverse')
set(gca, 'FontSize', 16); 
camlight('headlight') 
lighting gouraud

%-----------------------------------
subplot(2,3,4)
model = 41;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX/1000,YY/1000,ZZ/1000,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:)/1000,Xs(2,:)/1000,Xs(3,:)/1000, 'rv', 'Markersize', 10);
plot3(Xg(1,:)/1000,Xg(2,:)/1000,Xg(3,:)/1000, 'g^', 'Markersize', 10);
plot3(X0(1)/1000,X0(2)/1000,X0(3)/1000, 'b*', 'Markersize', 10);
axis([0 4 0 4 0 2])
xlabel('X  [km]', 'Fontsize', 16); 
ylabel('Y [km]','Fontsize', 16); 
zlabel('Depth  [km]','Fontsize', 16); 
title('Flat reflector','Fontsize', 28)
title('Sphere (R = 1 km)','Fontsize', 28)
set(gca,'xtick',[0 2 ])
set(gca,'ytick',[0 2 4])
set(gca,'ztick',[0 1 2])
view(3); 
set(gca, 'ZDir', 'reverse')
set(gca, 'FontSize', 16); 
camlight('headlight') 
lighting gouraud

%-----------------------------------
subplot(2,3,5)
model = 51;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX/1000,YY/1000,ZZ/1000,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:)/1000,Xs(2,:)/1000,Xs(3,:)/1000, 'rv', 'Markersize', 10);
plot3(Xg(1,:)/1000,Xg(2,:)/1000,Xg(3,:)/1000, 'g^', 'Markersize', 10);
plot3(X0(1)/1000,X0(2)/1000,X0(3)/1000, 'b*', 'Markersize', 10);
axis([0 4 0 4 0 2])
xlabel('X  [km]', 'Fontsize', 16); 
ylabel('Y [km]','Fontsize', 16); 
zlabel('Depth  [km]','Fontsize', 16); 
title('Ellipsoid','Fontsize', 28)
set(gca,'xtick',[0 2 ])
set(gca,'ytick',[0 2 4])
set(gca,'ztick',[0 1 2])
view(3); 
set(gca, 'ZDir', 'reverse')
set(gca, 'FontSize', 16); 
camlight('headlight') 
lighting gouraud

%-----------------------------------
subplot(2,3,6)
model = 61;
Get_model_parameters; 
Get_model_acquisition_geometry; 
[XX, YY] = meshgrid(G.xx, G.yy);
[ZZ, ind] = Get_model_surface(XX,YY,model);
surf(XX/1000,YY/1000,ZZ/1000,'FaceColor',[1 0 0],'EdgeColor','none'); 
hold on
plot3(Xs(1,:)/1000,Xs(2,:)/1000,Xs(3,:)/1000, 'rv', 'Markersize', 10);
plot3(Xg(1,:)/1000,Xg(2,:)/1000,Xg(3,:)/1000, 'g^', 'Markersize', 10);
plot3(X0(1)/1000,X0(2)/1000,X0(3)/1000, 'b*', 'Markersize', 10);
axis([0 4 0 4 0 2])
xlabel('X  [km]', 'Fontsize', 16); 
ylabel('Y [km]','Fontsize', 16); 
zlabel('Depth  [km]','Fontsize', 16); 
title('Complex surface','Fontsize', 28)
set(gca,'xtick',[0 2 ])
set(gca,'ytick',[0 2 4])
set(gca,'ztick',[0 1 2])
view(3); 
set(gca, 'ZDir', 'reverse')
set(gca, 'FontSize', 16); 
camlight('headlight') 
lighting gouraud















