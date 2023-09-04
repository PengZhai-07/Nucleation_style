%%
clear
clc
close all
%%
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[0 0 1000 1300]);%左下角位置，宽高

t = tiledlayout(4,3);
t.TileSpacing = 'tight';
t.Padding = 'tight';

nexttile(1,[2,3])
%%
gamma = pi/4;  % empirical constant parameter about geometry
mu = 3.20e10;  % Pa 
sigma = 40e6;  % Pa
a = 0.015;
a_b = 0.2:0.05:0.9;
b = a./a_b;
r = 1;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
T = zeros(1,length(a_b));
nT = zeros(1,length(a_b));

L = [0.4, 0.5,0.6,0.8,1,1.2,1.6,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125,160,200,250,300]*10^-3;

H = 0;    % half-width
NS = zeros(length(b),length(L),length(H));
Cohesive = zeros(length(b),length(L),length(H));
res_constraint = zeros(length(b),length(L),length(H));
% 
diff_NS = zeros(length(b),length(L),length(H));
EEP = zeros(length(b),length(L),length(H));
MSR = zeros(length(b),length(L),length(H));

W = 5000;    % unit:m
m = 0; mm=0;

for i = 1:length(b)
        % define the simulation time for different a/b
    if a_b(i)<=0.26
        T(i) = 2000;
    elseif (0.29<=a_b(i)) && (a_b(i)<=0.36)
        T(i) = 1200;
    elseif (0.39<=a_b(i)) && (a_b(i)<=0.49)
        T(i) = 900;
    elseif (0.49<=a_b(i)) && (a_b(i)<=0.61)
        T(i) = 600;
    elseif (0.64<=a_b(i)) && (a_b(i)<=0.76)
        T(i) = 300;
    else
        T(i) = 150;
    end
    mu_D = mu;  % Pa
    for j = 1:length(L)   
        for k = 1:length(H)
            if a_b(i) > 0.3781
                syms y 
                exp = W/y*tanh(2*gamma*H(k)/W*y+atanh(mu_D/mu)) -...
                   2/pi*mu_D*L(j)*b(i)/sigma/(b(i)-a)^2;       % Rubin and Ampuero for a/b>0.5
   %                 exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -... %
                    % pi/4*mu_D(i)*L(j)/sigma/(b-a); %                exp =
                    % 1/y*tanh(2*H(k)*gamma/W*y+atanh(mu_D/mu)) -... %
                    % mu_D*L(j)/sigma/(a/a_b(i)-a)/W;    % without pi/4?
                y = double(vpasolve(exp,[0,1000000])) ;
            else
                y = W/(2*1.3774*mu_D*L(j)/b(i)/sigma);                        % Rubin and Ampuero for a/b<0.3781
            end
            Ru(i,j,k) = y;
            NS(i,j,k) = W/y;       % meter
            kk = mu/W*2/pi;     % for antiplane shear strain with constant slip
            C_1(i,j,k) = b(i)/a*(1-kk*L(j)/b(i)/sigma);     % with equation (17) and kk= G*neta/L
            Cohesive(i,j,k) = 9*pi/32*mu_D*r*L(j)/b(i)./sigma;  % for fault zone cases, see IDINI's derivation
            res_constraint(i,j,k) = min(NS(i,j,k), Cohesive(i,j,k)); 
            res = 32;
            if (y>=1.0) && (res_constraint(i,j,k) > 400/80*3) 
                m = m+1;
                P_1(m,:) = [log10(L(j)*1000), a_b(i), L(j),b(i), T(i)];     % resolution is enough
            elseif (y<1.0) && (res_constraint(i,j,k) > 400/80*3) 
                mm = mm + 1;
                P_2(mm,:) = [log10(L(j)*1000), a_b(i), L(j),b(i), T(i)];  
            end
            % get the possible values from simulations
            filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
            '0.25_',num2str(res),'_',num2str(T(i)),'_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
            '/Maximum_slip_velocity.out'];
            
            if exist(filename, "file") == 0
                disp(filename)

                res = 48; 
                filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                '0.25_',num2str(res),'_300_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                '/Maximum_slip_velocity.out'];
                if exist(filename, "file") == 0
                    disp(filename)
                    
                    filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                    '0.25_',num2str(res),'_150_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                     '/Maximum_slip_velocity.out'];
                     if exist(filename, "file") == 0
                         disp(filename)
                           res = 64; 
                            filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                            '0.25_',num2str(res),'_75_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                            '/Maximum_slip_velocity.out'];
                            if exist(filename, "file") == 0
                                disp(filename)
                                res = 80; 
                                filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                                '0.25_',num2str(res),'_75_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                                '/Maximum_slip_velocity.out'];
                                if exist(filename, "file") == 0
                                    disp(filename)
                                    continue
                                end
                            end
                     end
                end
            end

            fileID  = fopen(filename,'r');
            % extract EEP
            A = textscan(fileID, '%f', 1, 'HeaderLines',0);
            temp_MSR = A{1}(1);
%             if temp_MSR > 100
%                 temp_MSR = 0;
%             end
            MSR(i,j,k) = temp_MSR;
            fclose(fileID);

        end
    end
end

%%
hold on
%%
[Y,X] = meshgrid(a_b, log10(L*1000));
A = pcolor(X,Y,log10(MSR'));
shading interp
% % mycolormap = customcolormap([0 0.4 0.7 1], [1 1 1;0 0 1; 0 1 0; 1 0 0]);
% % colormap(flipud(mycolormap))
colormap("turbo")
% colormap(flipud(othercolor('RdYlBu_11b')))
% clim([min(min(log(MSR))),max(max(log(MSR)))])
c = colorbar;
ylabel(c, 'Max. Sliprate(m/s)')
c.Ticks = [-9:3:2];
c.TickLabels = ["10^{-9}","10^{-6}","10^{-3}","1"];
clim([-9,2])

%% transition zone
plot([log10(0.5), log10(300)],[0.4,0.4], "k--")
plot([log10(0.5), log10(300)],[0.35,0.35], "k--")

%% no data region
n_xxx = [log10(1.6),log10(1.6),log10(1.2), log10(1.0), log10(0.8), log10(0.8), log10(0.8),log10(0.6),log10(0.6), log10(0.5), log10(0.5), log10(0.5),log10(0.4), log10(0.4)];
n_yyy = [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.2];
fill(n_xxx, n_yyy, [0.8,0.8,0.8],"LineStyle","-")
text(log10(1),0.3, 'No data','Rotation',0,'FontSize',12)
%%
% different rupture style
%  aseismic slip
j = [11:12,16:18,19:22,24:28, 29:35, 36:43, 45:53, 55:64, 66:77, 78:91, 93:108, 110:128, 130:151];
scatter(P_2(j,1), P_2(j,2),'o', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
%  bilateral
i = [20:21,39:41,59:61,79:81,99:104,119:125,140:145,160:165,180:183,196:201,213:217,227:232,240:245,250:255,257:262];
scatter(P_1(i,1), P_1(i,2),'h','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
j = [1:10, 13:15, 23, 54, 65, 92, 109, 129];
scatter(P_2(j,1), P_2(j,2),'h','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
%  unilateral
i = [13:19,33:38,53:58,72:78, 93:98, 116:118, 136:139, 157:159, 177:179, 194:195, 209:212, 223:226, 236:239, 247, 249, 256];
scatter(P_1(i,1), P_1(i,2),'square','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
%  full and partial
i = [10:12,31:32,51:52, 92, 115, 134:135,154:156,173:176,189:193,205:208,219:222,233:235,246, 248];
scatter(P_1(i,1), P_1(i,2),'v','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
%   crack-like
i = [8:9,29:30,50,70:71,90:91,112:114, 132:133, 152:153,170:172, 186:188, 203:204,218];
scatter(P_1(i,1), P_1(i,2),'p','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
%   pulse-like
i = [1:7, 22:28, 42:49, 62:69, 82:89, 105:111,126:131, 146:151,166:169,184:185,202];
scatter(P_1(i,1), P_1(i,2),'d','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
box on
%%
v = [1, 4, 8, 16, 32];
[c,h]=contour(X(:,5:end),Y(:,5:end),Ru(5:end,:)',v);
clabel(c,h)
set(h,"color","k")

%%
% set(h,"color", [1,1,1;0,0,0;])
% set(findobj(gca,'Type','patch','UserData',1),'EdgeColor', 'w')
v = [0.5 , 2, 8, 12, 18];
[c,h]=contour(X(:,1:4),Y(:,1:4),Ru(1:4,:)',v);
clabel(c,h)
set(h,"color","black")

%% representative cases
scatter(log10(63), 0.55, 100, 'o', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
scatter(log10(25), 0.55,100,'h','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
scatter(log10(8), 0.55,100,'square','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
scatter(log10(4), 0.55,100,'v','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
scatter(log10(2.5), 0.55,100,'p','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
scatter(log10(1), 0.55,100,'d','filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
%% rupture style for a/b>=0.4
text(log10(200),0.45,"Aseismic slip",'Color','k','Rotation',40)
text(log10(60),0.45,"Symmetric-bilateral",'Rotation',40)
text(log10(18),0.45,["Unsymmetric-";"bilateral";"and unilateral"],'Rotation',40)
text(log10(10),0.45,"Full and partial",'Rotation',40)
text(log10(5),0.45,"Crack-like with aftershocks ",'Rotation',40)
text(log10(2),0.45,["Pulse-like with";  "aftershocks"],'Rotation',40)

% text(log10(50),0.375,["Transition zone of two nucleation equations"],'Rotation',0)
% rupture style for a/b<=0.35
% text(log10(200),0.25,"SSE(<0.1m/s) and Creep(<1e-8m/s)",'Rotation',40)
text(log10(160),0.23,"Symmetric-bilateral",'Rotation',45)
text(log10(40),0.23,["Unsymmetric-";"bilateral";"and unilateral"],'Rotation',45)
text(log10(16),0.23,"Full and partial",'Rotation',45)
text(log10(10.2),0.23,["Crack-like with"; "aftershocks"],'Rotation', 55)
text(log10(5),0.23,["Pulse-like with";  "aftershocks"],'Rotation',45)

%%
xticks([log10(L*1000)])
xticklabels([0.4, 0.5,0.6,0.8,1,1.2,1.6,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125,160,200,250,300]/5)
yticks([0.2:0.05:0.9])
yticklabels([0.2:0.05:0.9])
set(gca,'XDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
%set(gca,'YDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
xlabel('RD_{RS}')
ylabel('a/b')
title("(a)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
%%
set(gca,'TickDir', 'out')

% axis normal
% print(gcf, 'Maximum_sliprate.png','-dpng','-r600')
% saveas(gcf, 'Maximum_sliprate.png')

%% slip velocity distribution

%% aseismic slip
nexttile(7)
filename = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_32_600_0_0_1.0_0.0_4_0.55_0.063/sliprate.out';
if exist(filename, "file") == 0
         disp(filename)
end     
sliprate = readmatrix(filename, "FileType","text");
ll = 6000; rr = 9000;
nn = 101:701;
threshold = 1e-4;
% down sampling
sliprate_1 = sliprate(ll:rr, :);
n = 1;m = 1;
indx_1 = zeros(3,1);
indx_2 = zeros(3,1);
for i = 2: rr-ll
    if max(sliprate_1(i-1,:))> threshold &&  max(sliprate_1(i,:))< threshold    % end of coseismic phase
        indx_2(n) = i;
        n = n + 1;
    end
    if max(sliprate_1(i-1,:)) < threshold&&  max(sliprate_1(i,:)) > threshold   % start of coseismic phase
        indx_1(m) = i;
        m = m + 1;
    end
end
new_v = [sliprate_1(1:10:indx_1(1),nn);sliprate_1(indx_1(1):indx_2(1),nn);sliprate_1(indx_2(1):10:indx_1(2),nn);...
    sliprate_1(indx_1(2):indx_2(2),nn);sliprate_1(indx_2(2):10:indx_1(3),nn); sliprate_1(indx_1(3):indx_2(3),nn);sliprate_1(indx_2(3):10:end,nn)];
pcolor(log10(new_v))
colormap('turbo')
shading interp
clim([-9,2])
xlabel('Distance (km)')
xticks([101 501])
xticklabels([-2.5 2.5])
% set(gca,'XTickLabel', [])
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
set(gca,'TickDir', 'out')
hold on
plot([101, 101],[0,rr-ll], 'k--')
plot([501, 501],[0,rr-ll], 'k--')
scatter(301, (rr-ll)*9/10, 80, 'o', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(201,(rr-ll)*4/5,"Aseismic slip",'Color','w')
title("(b)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
clear sliprate
%% bilateral
nexttile(8)
filename = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_32_600_0_0_1.0_0.0_4_0.55_0.025/sliprate.out';
if exist(filename, "file") == 0
         disp(filename)
end     
sliprate = readmatrix(filename, "FileType","text");
ll = 5000; rr = 15500;
nn = 101:701;
threshold = 1e-0;
% down sampling
sliprate_1 = sliprate(ll:rr, :);
n = 1;m = 1;
indx_1 = zeros(3,1);
indx_2 = zeros(3,1);
for i = 2: rr-ll
    if max(sliprate_1(i-1,:))> threshold &&  max(sliprate_1(i,:))< threshold    % end of coseismic phase
        indx_2(n) = i;
        n = n + 1;
    end
    if max(sliprate_1(i-1,:)) < threshold&&  max(sliprate_1(i,:)) > threshold   % start of coseismic phase
        indx_1(m) = i;
        m = m + 1;
    end
end
new_v = [sliprate_1(1:10:indx_1(1),nn);sliprate_1(indx_1(1):indx_2(1),nn);sliprate_1(indx_2(1):10:indx_1(2),nn);...
    sliprate_1(indx_1(2):indx_2(2),nn);sliprate_1(indx_2(2):10:indx_1(3),nn); sliprate_1(indx_1(3):indx_2(3),nn);sliprate_1(indx_2(3):10:end,nn)];
pcolor(log10(new_v))
colormap('turbo')
shading interp
clim([-9,2])
xlabel('Distance (km)')
xticks([101 501])
xticklabels([-2.5 2.5])
% set(gca,'XTickLabel', [])
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])

set(gca,'TickDir', 'out')
hold on
plot([101, 101],[0,rr-ll], 'k--')
plot([501, 501],[0,rr-ll], 'k--')
scatter(301, (rr-ll)*9/10, 80, 'h', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(151,(rr-ll)*4/5,"Symmetric-bilateral",'Color','w')
title("(c)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
clear sliprate
%% unilateral
nexttile(9)
filename = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_32_600_0_0_1.0_0.0_4_0.55_0.008/sliprate.out';
if exist(filename, "file") == 0
         disp(filename)
end     
sliprate = readmatrix(filename, "FileType","text", "OutputType","double");
ll = 5500; rr = 12500;
nn = 101:701;
threshold = 1e-0;
% down sampling
sliprate_1 = sliprate(ll:rr, :);
n = 1;m = 1;
indx_1 = zeros(3,1);
indx_2 = zeros(3,1);
for i = 2: rr-ll
    if max(sliprate_1(i-1,:))> threshold &&  max(sliprate_1(i,:))< threshold    % end of coseismic phase
        indx_2(n) = i;
        n = n + 1;
    end
    if max(sliprate_1(i-1,:)) < threshold&&  max(sliprate_1(i,:)) > threshold   % start of coseismic phase
        indx_1(m) = i;
        m = m + 1;
    end
end
new_v = [sliprate_1(1:10:indx_1(1),nn);sliprate_1(indx_1(1):indx_2(1),nn);sliprate_1(indx_2(1):10:indx_1(2),nn);...
    sliprate_1(indx_1(2):indx_2(2),nn);sliprate_1(indx_2(2):10:indx_1(3),nn); sliprate_1(indx_1(3):indx_2(3),nn);sliprate_1(indx_2(3):10:end,nn)];
pcolor(real(log10((new_v))))
colormap('turbo')
shading interp
clim([-9,1])
xlabel('Distance (km)')
xticks([101 501])
xticklabels([-2.5 2.5])
% set(gca,'XTickLabel', [])
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
set(gca,'TickDir', 'out')
hold on
plot([101, 101],[0,rr-ll], 'k--')
plot([501, 501],[0,rr-ll], 'k--')
scatter(301, (rr-ll)*9/10, 80, 'square', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(151,(rr-ll)*4/5,["Unsymmetric-bilateral";"and unilateral"],'Color','w')
title("(d)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
clear sliprate
%% full and partial
nexttile(10)
filename = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_32_600_0_0_1.0_0.0_4_0.55_0.004/sliprate.out';
if exist(filename, "file") == 0
         disp(filename)
end     
sliprate = readmatrix(filename, "FileType","text", "OutputType","double");
ll = 4100; rr = 7500;
nn = 101:701;
threshold = 1e0;
% down sampling
sliprate_1 = sliprate(ll:rr, :);
n = 1;m = 1;
indx_1 = zeros(4,1);
indx_2 = zeros(4,1);
for i = 2: rr-ll
    if max(sliprate_1(i-1,:))> threshold &&  max(sliprate_1(i,:))< threshold    % end of coseismic phase
        indx_2(n) = i;
        n = n + 1;
    end
    if max(sliprate_1(i-1,:)) < threshold&&  max(sliprate_1(i,:)) > threshold   % start of coseismic phase
        indx_1(m) = i;
        m = m + 1;
    end
end
new_v = [sliprate_1(1:10:indx_1(1),nn);sliprate_1(indx_1(1):indx_2(1),nn);sliprate_1(indx_2(1):10:indx_1(2),nn);...
    sliprate_1(indx_1(2):indx_2(2),nn);sliprate_1(indx_2(2):10:indx_1(3),nn); sliprate_1(indx_1(3):indx_2(3),nn);sliprate_1(indx_2(3):10:indx_1(4),nn);...
    sliprate_1(indx_1(4):indx_2(4),nn);sliprate_1(indx_2(4):10:end,nn)];
pcolor(real(log10((new_v))))
colormap('turbo')
shading interp
clim([-9,2])
xlabel('Distance (km)')
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
xticks([101 501])
xticklabels([-2.5 2.5])
set(gca,'TickDir', 'out')
hold on
plot([101, 101],[0,rr-ll], 'k--')
plot([501, 501],[0,rr-ll], 'k--')
scatter(301, (rr-ll)*9/10, 80, 'v', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(201,(rr-ll)*4/5,"Full and Partial",'Color','w')
title("(e)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
clear sliprate
%% crack-like with aftershocks
nexttile(11)
filename = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_32_600_0_0_1.0_0.0_4_0.55_0.0025/sliprate.out';
if exist(filename, "file") == 0
         disp(filename)
end     
sliprate = readmatrix(filename, "FileType","text", "OutputType","double");
ll = 5200; rr = 12500;
nn= 101:701;
threshold = 1e0;
%% down sampling
sliprate_1 = sliprate(ll:rr, :);
n = 1;m = 1;
indx_1 = zeros(11,1);
indx_2 = zeros(11,1);
for i = 2: rr-ll
    if max(sliprate_1(i-1,:))> threshold &&  max(sliprate_1(i,:))< threshold    % end of coseismic phase
        indx_2(n) = i;
        n = n + 1;
    end
    if max(sliprate_1(i-1,:)) < threshold&&  max(sliprate_1(i,:)) > threshold   % start of coseismic phase
        indx_1(m) = i;
        m = m + 1;
    end
end
new_v = [sliprate_1(1:10:indx_1(1),nn);sliprate_1(indx_1(1):indx_2(1),nn);sliprate_1(indx_2(1):10:indx_1(2),nn);...
    sliprate_1(indx_1(2):indx_2(2),nn);sliprate_1(indx_2(2):10:indx_1(3),nn); sliprate_1(indx_1(3):indx_2(3),nn);sliprate_1(indx_2(3):10:indx_1(4),nn);...
    sliprate_1(indx_1(4):indx_2(4),nn);sliprate_1(indx_2(4):10:indx_1(5),nn); sliprate_1(indx_1(5):indx_2(6),nn);sliprate_1(indx_2(6):10:indx_1(7),nn);...
    sliprate_1(indx_1(7):indx_2(7),nn);sliprate_1(indx_2(7):10:indx_1(8),nn); sliprate_1(indx_1(8):indx_2(8),nn);sliprate_1(indx_2(8):10:indx_1(9),nn);...
    sliprate_1(indx_1(9):indx_2(9),nn);sliprate_1(indx_2(9):10:indx_1(10),nn); sliprate_1(indx_1(10):indx_2(10),nn);sliprate_1(indx_2(10):10:indx_1(11),nn);...
    sliprate_1(indx_1(11):indx_2(11),nn);sliprate_1(indx_2(11):10:end,nn)];
pcolor(real(log10((new_v))))
colormap('turbo')
shading interp
clim([-9,2])
xlabel('Distance (km)')
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
xticks([101 501])
xticklabels([-2.5 2.5])
set(gca,'TickDir', 'out')
hold on
plot([101, 101],[0,rr-ll], 'k--')
plot([501, 501],[0,rr-ll], 'k--')
scatter(301, (rr-ll)*9/10, 80, 'p', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(151,(rr-ll)*4/5,"Crack-like with aftershocks",'Color','w')
title("(f)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
clear sliprate
%% pulse-like with afteshocks
nexttile(12)
filename = '/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.25_48_300_0_0_1.0_0.0_4_0.55_0.001/sliprate.out';
if exist(filename, "file") == 0
         disp(filename)
end     
sliprate = readmatrix(filename, "FileType","text", "OutputType","double");
ll = 9500; rr = 18000;
nn = 151:1051;
threshold = 1e-1;
% down sampling
sliprate_1 = sliprate(ll:rr, :);
n = 1;m = 1;
indx_1 = zeros(19,1);
indx_2 = zeros(19,1);
for i = 2: rr-ll
    if sliprate_1(i-1,601)> threshold &&  sliprate_1(i,601)< threshold    % end of coseismic phase
        indx_2(n) = i;
        n = n + 1;
    end
    if sliprate_1(i-1,601) < threshold &&  sliprate_1(i,601) > threshold   % start of coseismic phase
        indx_1(m) = i;
        m = m + 1;
    end
end
new_v = [sliprate_1(1:10:indx_1(1),nn);sliprate_1(indx_1(1):indx_2(1),nn);sliprate_1(indx_2(1):10:indx_1(2),nn);...
    sliprate_1(indx_1(2):indx_2(2),nn);sliprate_1(indx_2(2):10:indx_1(3),nn); sliprate_1(indx_1(3):indx_2(3),nn);sliprate_1(indx_2(3):10:indx_1(4),nn);...
    sliprate_1(indx_1(4):indx_2(4),nn);sliprate_1(indx_2(4):10:indx_1(5),nn); sliprate_1(indx_1(5):indx_2(6),nn);sliprate_1(indx_2(6):10:indx_1(7),nn);...
    sliprate_1(indx_1(7):indx_2(7),nn);sliprate_1(indx_2(7):10:indx_1(8),nn); sliprate_1(indx_1(8):indx_2(8),nn);sliprate_1(indx_2(8):10:indx_1(9),nn);...
    sliprate_1(indx_1(9):indx_2(9),nn);sliprate_1(indx_2(9):10:indx_1(10),nn);sliprate_1(indx_1(10):indx_2(10),nn);sliprate_1(indx_2(10):10:indx_1(11),nn);...
    sliprate_1(indx_1(11):indx_2(11),nn);sliprate_1(indx_2(11):10:indx_1(12),nn); sliprate_1(indx_1(12):indx_2(12),nn);sliprate_1(indx_2(12):10:indx_1(13),nn); ...
    sliprate_1(indx_1(13):indx_2(13),nn);sliprate_1(indx_2(13):10:indx_1(14),nn); sliprate_1(indx_1(14):indx_2(14),nn);sliprate_1(indx_2(14):10:indx_1(15),nn); ...
    sliprate_1(indx_1(15):indx_2(15),nn);sliprate_1(indx_2(15):10:indx_1(16),nn); sliprate_1(indx_1(16):indx_2(16),nn);sliprate_1(indx_2(16):10:indx_1(17),nn); ...
    sliprate_1(indx_1(17):indx_2(17),nn);sliprate_1(indx_2(17):10:indx_1(18),nn); sliprate_1(indx_1(18):indx_2(18),nn);sliprate_1(indx_2(18):10:indx_1(19),nn); ...
    sliprate_1(indx_1(19):10:end,nn)];
pcolor(real(log10((new_v))))
colormap('turbo')
shading interp
clim([-9,2])
xlabel('Distance (km)')
ylabel('Time Steps(Integer)')
set(gca,'YTickLabel', [])
xticks([151 751])
xticklabels([-2.5 2.5])
set(gca,'TickDir', 'out')
hold on
plot([151, 151],[0,rr-ll], 'k--')
plot([751, 751],[0,rr-ll], 'k--')
scatter(451, (rr-ll)*9/10, 80, 'd', 'filled','MarkerFaceColor','w', 'MarkerEdgeColor', 'k')
text(251,(rr-ll)*4/5,"Pulse-like with aftershocks",'Color','w')
title("(g)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
clear sliprate
%%
exportgraphics(gcf,'maximum_sliprate.png','Resolution',200)
% print -djpeg maximum_sliprate -r600
