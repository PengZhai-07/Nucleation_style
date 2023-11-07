clc;clear;
close all
%% zoom in of pulse-like characteristics
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[0 0 1000 900]);%左下角位置，宽高

t = tiledlayout(3,2);
t.TileSpacing = 'tight';
t.Padding = 'tight';

%% 
ax(1) = nexttile(1);
STF = 461;
dir = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_32_600_0_0_1.0_0.0_4_0.55_0.0025/';
if exist(dir, "dir") == 0
         disp(dir)
end     
sliprate = readmatrix(strcat(dir,"sliprate.out"), "FileType","text", "OutputType","double");
t_v = readmatrix(strcat(dir,"time_velocity.out"), "FileType","text", "OutputType","double");
t = t_v(:,1);
ll = 5280; rr = 5500;
pcolor(log10(abs(sliprate(ll:rr,101:701))));
colormap(ax(1),'hot')
shading interp
clim([-1,1])
% c = colorbar;
% ylabel(c, 'Slip Velocity (m/s)')
% c.Ticks = [-1:1:1];
% c.TickLabels = ["10^{-1}","1","10"];
xlabel('Distance (km)')
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
xticks([61 541])
xticklabels([-3 3])
set(gca,'TickDir', 'out')
hold on
plot([STF, STF],[0,rr-ll], 'b--')      % STF
scatter(101, (rr-ll)*9/10, 80, 'p', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(121,(rr-ll)*9/10,"Crack-like",'Color','w')
title("(a)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';

ax(3) = nexttile(3);
n_sliprate = abs(sliprate(ll:rr,101:701));
[a, b] = size(n_sliprate);
n_t = t(ll*10:10:rr*10)-t(ll*10);
n_t_diff = zeros(a-1);
slip = zeros(a-1, b);
cul_slip = zeros(a, b);
for i = 1:a-1
    n_t_diff(i) = n_t(i+1)-n_t(i);
    slip(i,:) = n_t_diff(i) .* (n_sliprate(i,:) + n_sliprate(i+1,:))./2;
    cul_slip(i+1,:) = cul_slip(i,:) + slip(i,:);
end
smooth_cul_slip_1 = smoothdata(cul_slip,2,"gaussian",50);
smooth_cul_slip = smoothdata(smooth_cul_slip_1, 1,"gaussian",50);
m_slip = max(max(cul_slip));
% extract the values of contours
x_interval = 0.005;
y_interval = 1;
v = x_interval:x_interval:m_slip;
YY = ones(length(n_t),1)*[1:b];      % slip
XX = n_t*ones(1,b);
[c,h]=contour(XX,YY,smooth_cul_slip,v);
CC = getContourLineCoordinates(c);
xq = 1:y_interval:601;
n = max(CC.Group);
m = max(CC.Level);
T = max(CC.X);
vq = zeros(length(xq),n);
for i = 1:n
    rows = (CC.Group==i);
    vq(:,i) = interp1(table2array(CC(rows,["Y"])),table2array(CC(rows,["X"])),xq,"linear");
end
xx = ones(length(xq),1)*[1:n].*x_interval;      % slip
yy = xq'.*ones(1,n);      % distance
vq1 = fillmissing(vq,  'movmean', 10, 1);
vq = fillmissing(vq1,  'movmean', 10, 2);
% vq(isnan(vq)) = T+1;
pcolor(yy,xx,vq);
shading interp
colormap(ax(3),flipud(turbo))
clim([0,T])
% c = colorbar;
% ylabel(c, 'Rupture time (s)')
% c.Ticks = [0:1:3];
% c.TickLabels = ["0","1","2","3"];
hold on
[c,h] = contour(yy,xx,vq,[0.2:0.2:T],  "LineWidth", 1.5);
set(h,"color","w")
xlabel('Distance (km)')
ylabel('Slip (m)')
xticks([61 541])
xticklabels([-3 3])
set(gca,'TickDir', 'out')
title("(c)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
axis([0 601 0 1.5])

nexttile(5)
plot(n_t, sliprate(ll:rr,STF), "LineWidth",1.5)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
axis([0 3.5 0 12])
title("(e)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';

clear sliprate

%%
ax(2) = nexttile(2);
STF = 681;
dir = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_48_300_0_0_1.0_0.0_4_0.55_0.001/';
if exist(dir, "dir") == 0
         disp(dir)
end     
sliprate = readmatrix(strcat(dir,"sliprate.out"), "FileType","text", "OutputType","double");
t_v = readmatrix(strcat(dir,"time_velocity.out"), "FileType","text", "OutputType","double");
t = t_v(:,1);
ll = 9570; rr = 9900;
pcolor(log10(abs(sliprate(ll:rr,151:1051))))
colormap(ax(2),'hot')
shading interp
clim([-1,1])
c = colorbar;
ylabel(c, 'Slip velocity (m/s)')
c.Ticks = [-1:1:1];
c.TickLabels = ["10^{-1}","1","10"];
xlabel('Distance (km)')
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
xticks([91 811])
xticklabels([-3 3])
set(gca,'TickDir', 'out')
hold on
plot([STF, STF],[0,rr-ll], 'b--')    % points for the STF
scatter(151, (rr-ll)*9/10, 80, 'd', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(181,(rr-ll)*9/10,["Combination of"; "pulse and crack"],'Color','w')
title("(b)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';

ax(4) = nexttile(4);
n_sliprate = abs(sliprate(ll:rr,151:1051));
[a, b] = size(n_sliprate);
n_t = t(ll*10:10:rr*10)-t(ll*10);
n_t_diff = zeros(a-1);
slip = zeros(a-1, b);
cul_slip = zeros(a, b);
for i = 1:a-1
    n_t_diff(i) = n_t(i+1)-n_t(i);
    slip(i,:) = n_t_diff(i) .* (n_sliprate(i,:) + n_sliprate(i+1,:))./2;
    cul_slip(i+1,:) = cul_slip(i,:) + slip(i,:);
end
smooth_cul_slip_1 = smoothdata(cul_slip,2,"gaussian",50);
smooth_cul_slip = smoothdata(smooth_cul_slip_1,1,"gaussian",50);
m_slip = max(max(cul_slip));
% extract the values of contours
x_interval = 0.005;
y_interval = 1;
v = x_interval:x_interval:m_slip;
YY = ones(length(n_t),1)*[1:b];      % slip
XX = n_t*ones(1,b);
[c,h]=contour(XX,YY,smooth_cul_slip,v);
CC = getContourLineCoordinates(c);
xq = 1:y_interval:901;
n = max(CC.Group);
m = max(CC.Level);
T = max(CC.X);
vq = zeros(length(xq),n);
for i = 1:n
    rows = (CC.Group==i);
    vq(:,i) = interp1(table2array(CC(rows,["Y"])),table2array(CC(rows,["X"])),xq,"linear");
end
xx = ones(length(xq),1)*[1:n].*x_interval;      % slip
yy = xq'.*ones(1,n);
% vq1 = inpaint_nans(vq, 1);
vq1 = fillmissing(vq,  'movmean', 10, 1);
vq = fillmissing(vq1,  'movmean', 10, 2);
% vq(isnan(vq)) = T+1;
pcolor(yy,xx,vq);
shading interp
colormap(ax(4),flipud(turbo))
clim([0,T])
c = colorbar;
ylabel(c, 'Rupture time (s)')
c.Ticks = [0:1:3];
c.TickLabels = ["0","1","2","3"];
hold on
[c,h] = contour(yy,xx,vq,[0.2:0.2:T], "LineWidth", 1.5);
set(h,"color","w")
xlabel('Distance (km)')
ylabel('Slip (m)')
xticks([91 811])
xticklabels([-3 3])
set(gca,'TickDir', 'out')
title("(d)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
axis([0 901 0 1.2])

nexttile(6)
plot(n_t, sliprate(ll:rr,STF), "LineWidth",1.5)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
axis([0 3.5 0 12])
title("(f)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';

clear sliprate


%%
exportgraphics(gcf,'Comparison_crack_pulse.png','Resolution',600)
