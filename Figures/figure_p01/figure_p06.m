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
cdpmin = 2400;
cdpmax = 2800;
cdpcmp = 2500;
tt = (0:2500)*0.002; 
cc = (0:333)*0.02; 


%% Upload original prestack data

data = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_cdp_2500_proc.su');

header = get_segy_headers(data);

%% Plot CMP gather
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

subplot(1,5,[1 2])

ind = header.cdp == cdpcmp; 
CMP = data.traces(:,ind);
h = abs(header.offset(ind)/2000);
% imagesc(h,tt,CMP);
% colormap('gray')
% set(gca, 'YDir', 'reverse')
% shading interp
% hold on
% caxis([-4 4])

for i = 1:length(h)
    f = CMP(:,i);
    f1 = 0.2*f; 
    f2 = 0.4*f; 
    f3 = 0.6*f; 
    f4 = 0.8*f; 
    p = f < 0; 
    f1(p) = f(p);
    f2(p) = f(p);
    f3(p) = f(p);
    f4(p) = f(p);
    hold on
    plot(h(i)+0.002*f1, tt, 'black','LineWidth', 1);
    plot(h(i)+0.002*f2, tt, 'black','LineWidth', 1);
    plot(h(i)+0.002*f3, tt, 'black','LineWidth', 1);
    plot(h(i)+0.002*f4, tt, 'black','LineWidth', 1);
    plot(h(i)+0.002*f, tt, 'black','LineWidth', 1 );

end

if (1)
    t0 = 1.40; 
    hold on 
    vNMO = 1.480; 
    hh = 0:0.001:max(h); 
    t = sqrt(t0^2 + 4*hh.^2/vNMO^2); 
    plot(hh, t, '--r', 'LineWidth', 3); 
    
    t0 = 2.18; 
    vNMO = 1.600; 
    t = sqrt(t0^2 + 4*hh.^2/vNMO^2); 
    plot(hh, t, '--r', 'LineWidth', 3); 
    plot(tt*0, tt, '--b', 'LineWidth', 3); 
    
    
    t0 = 2.65; 
    vNMO = 1.480; 
    t = sqrt(t0^2 + 4*hh.^2/vNMO^2); 
    plot(hh, t, '--r', 'LineWidth', 3); 
end
plot(tt(1:10:end)*0, tt(1:10:end), '--b', 'LineWidth', 3); 
set(gca, 'YDir', 'reverse');
xlabel('Half-offset, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 0.5 1.0 1.5], 'FontSize',30) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30 )
title('CMP gather','FontSize',30)
axis([0 1.9 1 4])

 
%% Upload stacked data

data_cmp = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs.DE.Fstack_proc.su');


headers_cmp = get_segy_headers(data_cmp); 
[cdp, ind] = sort(headers_cmp.cdp);
traces_cmp = data_cmp.traces(:,ind); 


%% Plot stacked trace

ind = cdp == cdpcmp;

subplot(1,5,3)
plot(SME(traces_cmp(:,ind)/10,1), tt, 'black','LineWidth', 1)
axis([-1 1 1 4])
set(gca, 'YDir', 'reverse'); 
xlabel('Amplitude','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[-1 0 1], 'FontSize',30) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30)
title('Stacked trace','FontSize',30)


%% Upload stacked section
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind3 = cdp == cdpcmp;
ind = (ind1+ind2) == 2; 

cc = headers_cmp.sx(ind);
ccx = abs(headers_cmp.sx(ind3)-cc(1))/1000;
cc = abs(cc - cc(1))/1000;

subplot(1,5,[4 5])
imagesc(cc, tt, traces_cmp(:,ind)/5)
colormap('gray')
title('CMP stack','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2 4 6], 'FontSize',30) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([-1 1])
c = colorbar('Ticks', [-1, 0, 1]);
c.Label.String = 'Amplitude';
hold on
plot(linspace(ccx,ccx,251),tt(1:10:end), '--b', 'LineWidth', 3);

%% Plot stacked section in color 

colormap(makeColorMap([0 0 1],[1 1 1],[1 0 0],100));

%% Upload coherence stack data

data_coher = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_coher.su');


headers_coher = get_segy_headers(data_coher); 
[cdp, ind] = sort(headers_coher.cdp);
traces_coher = data_coher.traces(:,ind); 

%% Plot coherency
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind3 = cdp == cdpcmp;
ind = (ind1+ind2) == 2; 

subplot(1,5,[4 5])
imagesc(cc, tt, traces_coher(1:2501,ind))
caxis([-5 5])
%colormap('gray')
colormap('hot')
title('Coherency','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2 4 6], 'FontSize',30 ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30 )
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
hold on
caxis([0.0 1.0])
c = colorbar('Ticks', [0, 0.5, 1]);
c.Label.String = 'Coherency';
plot(linspace(ccx,ccx,251),tt(1:10:end), '--b', 'LineWidth', 3);

%% Upload velocity

data_vnmo = read_su_file('/scratch/local1/ivan/CRS_data/TGS/tgs_vnmo.su');


headers_vnmo = get_segy_headers(data_vnmo); 
[cdp, ind] = sort(headers_vnmo.cdp);
traces_vnmo = data_vnmo.traces(:,ind); 

%% Plot velocity
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind3 = cdp == cdpcmp;
ind = (ind1+ind2) == 2; 

subplot(1,5,[4 5])
imagesc(cc, tt, traces_vnmo(1:2501,ind)/1000)
caxis([-5 5])
colormap('jet')
title('Velocity','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2 4 6], 'FontSize',30) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30)
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
hold on
caxis([1 2])
c = colorbar('Ticks', [1, 1.5, 2]);
c.Label.String = 'Velocity, km/s';
plot(linspace(ccx,ccx,251),tt(1:10:end), '--b', 'LineWidth', 3);

%% Part 2

figure(2)
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

ax1=subplot(1,3,1)
imagesc(cc, tt, traces_cmp(:,ind)/5)
colormap(ax1,'gray')
title('CMP stack','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2 4 6], 'FontSize',30) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30)
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
caxis([-1 1])
c = colorbar('southoutside','Ticks', [-1, 0, 1]);
c.Label.String = 'Amplitude';


ax2=subplot(1,3,2)
imagesc(cc, tt, traces_coher(1:2501,ind))
colormap(ax2,'hot')
title('Coherency','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2 4 6], 'FontSize',30 ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30)
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
hold on
caxis([0.0 1.0])
c = colorbar('southoutside','Ticks', [0, 0.5, 1]);
c.Label.String = 'Coherency';

ax3=subplot(1,3,3)
imagesc(cc, tt, traces_vnmo(1:2501,ind)/1000)
colormap(ax3,'jet')
title('Velocity','FontSize',30)
xlabel('Distance, km','FontSize',30)
ylabel('Time, s','FontSize',30)
set(gca,'xtick',[0 2 4 6], 'FontSize',30) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',30)
axis([min(cc) max(cc) 1.0 4.0])
set(gca, 'YDir', 'reverse')
shading interp
hold on
caxis([1 2])
c = colorbar('southoutside','Ticks', [1, 1.5, 2]);
c.Label.String = 'Velocity, km/s';


%% Part 3 Adaptive filtering

