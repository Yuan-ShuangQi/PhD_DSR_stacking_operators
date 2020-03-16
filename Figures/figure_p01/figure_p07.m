%% Make figure with idea of CMP stack
%
% *Author*: Abakumov Ivan
%
% *Publication date*: 14th January 2017

%% Add MLIB

clear all; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;

%% Data settings
cdpmin = 2351;
cdpmax = 4500;
cdpcmp = 2500;
tt = (0:2500)*0.002; 
cc = (0:333)*0.02; 

%% Upload stacked data

data_cmp = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs.DE.Fstack_proc.su');

headers_cmp = get_segy_headers(data_cmp); 
[cdp, ind] = sort(headers_cmp.cdp);
traces_cmp = data_cmp.traces(:,ind); 

%% Plot CMP gather
figure(1)
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind3 = cdp == cdpcmp;
ind = (ind1+ind2) == 2; 

cc = headers_cmp.sx(ind);
ccx = abs(headers_cmp.sx(ind3)-cc(1))/1000;
cc = abs(cc - cc(1))/1000;

pcolor(cc, tt, traces_cmp(:,ind)/5)
colormap('gray')
title('CMP stack','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 5 10 15 20 25], 'FontSize',30 ) 
%set(gca,'xtick',[0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25], 'FontSize',24, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
%axis([min(cc) 12.5 1.0 4.0])
%axis([12.5 max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([-1 1])
colormap('gray')
%colormap(makeColorMap([0 0 1],[1 1 1],[1 0 0],100));

%% Add close-up 1 
hold on
plot(linspace(cc(1),cc(1),126),tt(501:10:1751), '--r', 'LineWidth', 3);
plot(linspace(cc(801),cc(801),126),tt(501:10:1751), '--r', 'LineWidth', 3);
plot(linspace(cc(1),cc(801),126),linspace(tt(501),tt(501),126), '--r', 'LineWidth', 3);
plot(linspace(cc(1),cc(801),126),linspace(tt(1751),tt(1751),126), '--r', 'LineWidth', 3);

%% Add close-up 2 
hold on
plot(linspace(cc(1450),cc(1450),101),tt(1501:5:2001), '--r', 'LineWidth', 3);
plot(linspace(cc(2050),cc(2050),101),tt(1501:5:2001), '--r', 'LineWidth', 3);
plot(linspace(cc(1450),cc(2050),126),linspace(tt(1501),tt(1501),126), '--r', 'LineWidth', 3);
plot(linspace(cc(1450),cc(2050),126),linspace(tt(2001),tt(2001),126), '--r', 'LineWidth', 3);

%% Add description of multiple
hold on
text(3,3.6,'Multiple','FontSize',30,'Color','red')
quiver(0.5,3.0,1,-0.35,0,'Color', 'red','LineWidth',4)
quiver(5,3.5,1,-0.35,0,'Color', 'red','LineWidth',4)
quiver(10,3.6,1,-0.35,0,'Color', 'red','LineWidth',4)
quiver(15,3.85,1,-0.35,0,'Color', 'red','LineWidth',4)
quiver(20,3.9,1,-0.35,0,'Color', 'red','LineWidth',4)
quiver(25,3.95,1,-0.35,0,'Color', 'red','LineWidth',4)

%% Add description of multiples
hold on
text(18,3.0,'Multiples','FontSize',30,'Color','red')
quiver(15,3.05,1,0.35,0,'Color', 'red','LineWidth',4)
quiver(20,3.1,1,0.35,0,'Color', 'red','LineWidth',4)
quiver(25,3.15,1,0.35,0,'Color', 'red','LineWidth',4)

%%
data_bepic = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_pic_bef_proc.su');
headers_bepic = get_segy_headers(data_bepic); 
[cdp_pic, ind] = sort(headers_bepic.cdp);
traces_bepic = data_bepic.traces(:,ind); 
cc_pic = cc(1:801);

data_afpic = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_pic_after_proc.su');
headers_afpic = get_segy_headers(data_afpic); 
[~, ind] = sort(headers_afpic.cdp);
traces_afpic = data_afpic.traces(:,ind); 

%% Plot
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 
pcolor(cc_pic, tt(501:1751), traces_afpic/5)
%pcolor(cc_pic, tt(501:1751), traces_bepic/5)
title('CMP stack','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2.5 5 7.5 10], 'FontSize',30) 
set(gca,'ytick',[1 2 3 4], 'FontSize',30)
axis([cc(1) cc(801) 1.0 3.5])
set(gca, 'YDir', 'reverse')
shading interp
caxis([-1 1])
colormap('gray')

%% Add multiple
hold on
text(3.5,3.33,'Multiple','FontSize',30,'Color','red')
quiver(0.5,2.9,1,-0.25,0,'Color', 'red','LineWidth',4)
quiver(5,3.4,1,-0.25,0,'Color', 'red','LineWidth',4)


%% MS with velocity
figure(3)
data_afvel = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_vel_after_proc.su');
headers_afvel = get_segy_headers(data_afvel); 
[cdp_vel, ind] = sort(headers_afvel.cdp);
traces_afvel = data_afvel.traces(:,ind); 
traces_afvel(1:101,1:240) = 0;
traces_afvel(1:226,241:962) = 0;
traces_afvel(1:301,963:end) = 0;

ind1 = cdp_vel >= cdpmin; 
ind2 = cdp_vel <= cdpmax; 
ind = (ind1+ind2) == 2; 
traces_afvel = traces_afvel(:,ind); 

cc_vel = cc;

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 
pcolor(cc, tt(501:2152), traces_afvel/5)
%pcolor(cc_pic, tt(501:1751), traces_bepic/5)
title('CMP stack','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 5 10 15 20 25], 'FontSize',30) 
set(gca,'ytick',[1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([-1 1])
%colormap(makeColorMap([0 0 1],[1 1 1],[1 0 0],100));
colormap('gray')


%% MS with coherence
figure(4)


data_afcoh = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_coher_after_proc.su');
headers_afcoh = get_segy_headers(data_afcoh); 
[cdp_coh, ind] = sort(headers_afcoh.cdp);
traces_afcoh = data_afcoh.traces(:,ind); 

ind1 = cdp_coh >= cdpmin; 
ind2 = cdp_coh <= cdpmax; 
ind = (ind1+ind2) == 2; 
traces_afcoh = traces_afcoh(501:2001,ind); 
traces_afcoh(1:101,1:240) = 0;
traces_afcoh(1:226,241:962) = 0;
traces_afcoh(1:301,963:end) = 0;
cdp_coh = cdp_coh(ind);

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 
pcolor(cc, tt(501:2001), traces_afcoh/5)
%pcolor(cc_pic, tt(501:1751), traces_bepic/5)
title('CMP stack','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 5 10 15 20 25], 'FontSize',30 ) 
set(gca,'ytick',[1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([-1 1])
%colormap(makeColorMap([0 0 1],[1 1 1],[1 0 0],100));
colormap('gray')

%% Attenuation using velocity (FIG 98)

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

subplot(1,2,1)
imagesc(cc(1450:2050), tt(1501:2001), traces_cmp(1501:2001,1450:2050)/5)
title(' Before','FontName', 'Arial','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[20 22.5 25], 'FontSize',30 ) 
set(gca,'ytick',[3 3.5 4], 'FontSize',30 )
axis([cc(1450) cc(2050) 3.0 4.0])
set(gca, 'YDir', 'reverse')
set(gca,'TickLabelInterpreter','latex')
shading interp
caxis([-1 1])
colormap('gray')

subplot(1,2,2)
imagesc(cc(1450:2050), tt(1501:2001), traces_afvel(1001:1501,1450:2050)/5)
title('After','FontName', 'Arial','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[20 22.5 25], 'FontSize',30 ) 
set(gca,'ytick',[3 3.5 4], 'FontSize',30 )
axis([cc(1450) cc(2050) 3.0 4.0])
set(gca, 'YDir', 'reverse')
set(gca,'TickLabelInterpreter','latex')
shading interp
caxis([-1 1])
colormap('gray')

%% Attenuation using coherency (FIG 88)

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

subplot(1,2,1)
imagesc(cc(1450:2050), tt(1501:2001), traces_cmp(1501:2001,1450:2050)/5)
title('Before','FontName', 'Arial','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[20 22.5 25], 'FontSize',30 ) 
set(gca,'ytick',[3 3.5 4], 'FontSize',30 )
axis([cc(1450) cc(2050) 3.0 4.0])
set(gca, 'YDir', 'reverse')
set(gca,'TickLabelInterpreter','latex')
shading interp
caxis([-1 1])
colormap('gray')

subplot(1,2,2)
imagesc(cc(1450:2050), tt(1501:2001), traces_afcoh(1001:1501,1450:2050)/5)
title('After','FontName', 'Arial','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[20 22.5 25], 'FontSize',30 ) 
set(gca,'ytick',[3 3.5 4], 'FontSize',30 )
axis([cc(1450) cc(2050) 3.0 4.0])
set(gca, 'YDir', 'reverse')
set(gca,'TickLabelInterpreter','latex')
shading interp
caxis([-1 1])
colormap('gray')

%% Add arrows
hold on
quiver(20.85, 3.65, -0.25,0.1,0,'Color', 'red','LineWidth',4)
quiver(21.45, 3.60, -0.25,0.1,0,'Color', 'red','LineWidth',4)
quiver(23.75, 3.47, -0.25,0.1,0,'Color', 'red','LineWidth',4)

%% Plot muted coherency section 

data_coh = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_coher_mute.su');
headers_coh = get_segy_headers(data_coh); 
[cdp_mcoh, ind] = sort(headers_coh.cdp);
traces_coh = data_coh.traces(501:2001,ind); 
traces_coh = traces_coh(:,1:2150); 

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

pcolor(cc, tt(501:2001), traces_coh)
colormap('hot')
title('Coherency','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 5 10 15 20 25], 'FontSize',30 ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([0.0 1.0])
c = colorbar('Ticks', [0, 0.5, 1]);
c.Label.String = 'Coherency';

%% Plot muted velocity section 

data_vnmo = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_vnmo_upmute.su');
headers_vnmo = get_segy_headers(data_vnmo); 
[cdp_mvnmo, ind] = sort(headers_vnmo.cdp);
traces_vnmo = data_vnmo.traces(501:2001,ind); 

fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

pcolor(cc, tt(501:2001), traces_vnmo/1000)

%%
colormap('jet')
title('Velocity','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 5 10 15 20 25], 'FontSize',30 ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([1 2])
c = colorbar('Ticks', [1, 1.5, 2]);
c.Label.String = 'Velocity';








