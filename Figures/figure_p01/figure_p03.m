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
cdpmin = 181903;
cdpmax = 182236;
cdpcmp = 182150;
tt = (0:624)*0.008; 
cc = (0:333)*0.02; 


%% Upload original prestack data

data = read_su_file('/home/zmaw/u250128/Desktop/CRS/SEG_C3_WA/3D_SEG_cdp_181903_182236.su');
header = get_segy_headers(data);

%% Plot CMP gather
fig = gcf;
fig.Color = [1.0 1.0 243/255];
fig.InvertHardcopy = 'off'; 

subplot(1,5,[1 2])

ind = header.cdp == cdpcmp; 
CMP = data.traces(:,ind);
h = header.offset(ind)/2000;
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
    plot(h(i)+0.015*f1, tt, 'black','LineWidth', 1);
    plot(h(i)+0.015*f2, tt, 'black','LineWidth', 1);
    plot(h(i)+0.015*f3, tt, 'black','LineWidth', 1);
    plot(h(i)+0.015*f4, tt, 'black','LineWidth', 1);
    plot(h(i)+0.015*f, tt, 'black','LineWidth', 1 );

end

if (1)
    t0 = 1.55; 
    hold on 
    vNMO = 2.100; 
    hh = 0:0.001:max(h); 
    t = sqrt(t0^2 + 4*hh.^2/vNMO^2); 
    plot(hh, t, '--r', 'LineWidth', 3); 
    vNMO = 1.500; 
    t = sqrt(t0^2 + 4*hh.^2/vNMO^2); 
    plot(hh, t, '--r', 'LineWidth', 3); 

    plot(tt*0, tt, '--b', 'LineWidth', 3); 
end
set(gca, 'YDir', 'reverse'); 
xlabel('Half-offset, km','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[0 0.5 1.0 1.5], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
title('CMP gather','FontSize',25, 'Color',[51/255 51/255 153/255])
axis([0 1.4 0.16 4.5])

 
%% Upload stacked data

data_cmp = read_su_file('/scratch/local1/ivan/CRS_data/SEG_C3_WA_3D_CRS/product_SEG_C3WA_ffid_proc_red.su_323/cmp/cmpsa.stack.su'); 
headers_cmp = get_segy_headers(data_cmp); 
[cdp, ind] = sort(headers_cmp.cdp);
traces_cmp = data_cmp.traces(:,ind); 


%% Plot stacked trace

ind = cdp == cdpcmp;

subplot(1,5,3)
plot(SME(traces_cmp(:,ind)/5,1), tt, 'black','LineWidth', 1)
axis([-1 1 0.16 4.5])
set(gca, 'YDir', 'reverse'); 
xlabel('Amplitude','FontSize',25,'Color',[51/255 51/255 153/255])
ylabel('Time, s','FontSize',25, 'Color',[51/255 51/255 153/255])
set(gca,'xtick',[-1 0 1], 'FontSize',20, 'XColor',[51/255 51/255 153/255] ) 
set(gca,'ytick',[0 1 2 3 4], 'FontSize',20,'YColor',[51/255 51/255 153/255] )
title('Stacked trace','FontSize',25, 'Color',[51/255 51/255 153/255])


%% Upload stacked section
ind1 = cdp >= cdpmin; 
ind2 = cdp <= cdpmax; 
ind = (ind1+ind2) == 2; 

subplot(1,5,[4 5])
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
ccx = cc(cdpcmp-cdpmin); 
hold on
plot(linspace(ccx,ccx,625),tt, '--b', 'LineWidth', 3);
