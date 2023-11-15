%%
clear
clc
close all
%%
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[0 0 1000 650]);%左下角位置，宽高
hold on
%%
vv = 0.25;
gamma = pi/4;  % empirical constant parameter about geometry
mu = 3.2e10;  % Pa 
sigma = 40e6;  % Pa
a = 0.015;
a_b = 0.2:0.05:0.9;
b = a./a_b;
r = 1;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
T = zeros(1,length(a_b));
nT = zeros(1,length(a_b));

L = [0.4,0.5,0.6,0.8,1.0,1.2,1.6,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125,160,200,250,300]*10^-3;

H = 0;    % half-width
NS = zeros(length(b),length(L),length(H));
Cohesive = zeros(length(b),length(L),length(H));
res_constraint = zeros(length(b),length(L),length(H));
% 
diff_NS = zeros(length(b),length(L),length(H));
EEP = zeros(length(b),length(L),length(H));

W = 5000;    % unit:m
m = 0;

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
%                 exp = W/y*tanh(2*gamma*H(k)/W*y+atanh(mu_D/mu)) -...
%                    2/pi*mu_D*L(j)*b(i)/sigma/(b(i)-a)^2/(1-vv);       % Rubin and Ampuero for a/b>0.5 and mode 2
                exp = W/y*tanh(2*gamma*H(k)/W*y+atanh(mu_D/mu)) -...
                   2/pi*mu_D*L(j)*b(i)/sigma/(b(i)-a)^2;       % Rubin and Ampuero for a/b>0.5 and mode 3
   %                 exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -... %
                    % pi/4*mu_D(i)*L(j)/sigma/(b-a); %                exp =
                    % 1/y*tanh(2*H(k)*gamma/W*y+atanh(mu_D/mu)) -... %
                    % mu_D*L(j)/sigma/(a/a_b(i)-a)/W;    % without pi/4?
                y = double(vpasolve(exp,[0,1000000])) ;
            else
                y = W/(2*1.3774*mu_D*L(j)/b(i)/sigma/(1-vv));                        % Rubin and Ampuero for a/b<0.3781
            end
            Ru(i,j,k) = y;
            NS(i,j,k) = W/y;       % meter
            kk = mu/W*2/pi;     % for antiplane shear strain with constant slip
            C_1(i,j,k) = b(i)/a*(1-kk*L(j)/b(i)/sigma);     % with equation (17) and kk= G*neta/L
            Cohesive(i,j,k) = 9*pi/32*mu_D*r*L(j)/b(i)./sigma;  % for fault zone cases, see IDINI's derivation
            res_constraint(i,j,k) = min(NS(i,j,k), Cohesive(i,j,k)); 
            if (y>=1.0) && (res_constraint(i,j,k) > 400/80*3) 
                m = m+1;
                P_1(m,:) = [log10(L(j)*1000), a_b(i), L(j),b(i), T(i)];     % resolution is enough
            end
            % get the possible values from simulations
            res = 32;
            filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
            '0.25_',num2str(res),'_',num2str(T(i)),'_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
            '/nucleation_stats.out'];
            
            if exist(filename, "file") == 0
                disp(filename)
                res = 48; 
                filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                '0.25_',num2str(res),'_300_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                '/nucleation_stats.out'];
                if exist(filename, "file") == 0
                    disp(filename)
                    filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                    '0.25_',num2str(res),'_150_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                     '/nucleation_stats.out'];
                     if exist(filename, "file") == 0
                         disp(filename)
                           res = 64; 
                            filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                            '0.25_',num2str(res),'_75_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                            '/nucleation_stats.out'];
                            if exist(filename, "file") == 0
                                disp(filename)
                                res = 80; 
                                filename = ['/home/pengzhai/Spear-normal-stress/plots/wholespace/phase_diagram_L_b/',...
                                '0.25_',num2str(res),'_75_0_0_1.0_0.0_4_',num2str(a_b(i)),'_',num2str(L(j)),...
                                '/nucleation_stats.out'];
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
            A = textscan(fileID, '%f', 4, 'HeaderLines',3);
            temp_EEP = A{1}(1);
            EEP(i,j,k) = temp_EEP;
            % different colors represent different value of EEP
%             if EEP(i,j,k) > 0.4    % constant weakening: fracture
%                 scatter(log10(1000*L(j)), a_b(i),60,'o','MarkerEdgeColor',[1,0,0], 'MarkerFaceColor',[1,0,0])
%             elseif   (EEP(i,j,k) <= 0.4) && (EEP(i,j,k) >0.2)    
%                 scatter(log10(1000*L(j)), a_b(i),60,'o','MarkerEdgeColor',[1,0.5,0], 'MarkerFaceColor',[1,0.5,0] )
%             elseif   (EEP(i,j,k) <= 0.2) && (EEP(i,j,k) >-0.2)    
%                 scatter(log10(1000*L(j)), a_b(i),60,'o','MarkerEdgeColor',[0,1,0],  'MarkerFaceColor',[0,1,0])
%             else  
%                 scatter(log10(1000*L(j)), a_b(i),60,'o','MarkerEdgeColor',[0,0,1],  'MarkerFaceColor',[0,0,1])
%             end
            
            fclose(fileID);

        end
    end
end

[Y,X] = meshgrid(a_b, log10(L*1000));
% A = pcolor(X,Y,EEP');
% set(A,'EdgeColor', 'None')
% shading interp

%%
mycolormap = customcolormap([0 0.28 0.32 1], [0 0 1; 1 1 1; 1 1 1; 1 0 0]);
clr = colormap(flipud(mycolormap));
% clim([min(min(EEP)), max(max(EEP))])
clim([-0.3, 0.7])
c = colorbar;
% c = colorbar('Ticks', linspace(-0.3, 0.6, 10), ...
%     'TickLabels', ["-0.3" "-0.2" "-0.1" "0.0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6"]);
ylabel(c, 'EEP'),
%%
EEP_T = EEP';
[i,j,v] = find(EEP_T);
for k = 1:length(v)
    n = EEP_T(i(k),j(k));      % EEP value     
    if EEP_T(i(k),j(k)) > 0.7
        n = 0.7;
    end
    if EEP_T(i(k),j(k)) < -0.3
        n = -0.299;
    end
    ratio = ceil((n + 0.3)*256); 
    scatter(X(i(k),j(k)),Y(i(k),j(k)), 80, 'O','filled','MarkerFaceColor',clr(ratio,:), 'MarkerEdgeColor', 'none')
end

grid off
hold on 

%%
v = [1,4];
[c,h]=contour(X(:,5:end),Y(:,5:end),Ru(5:end,:)',v);
clabel(c,h)
set(h,"color","black")
hold on
% v = [0.5, 2, 8, 16, 28];
v = [1,2, ];
[c,h]=contour(X(:,1:4),Y(:,1:4),Ru(1:4,:)',v);
clabel(c,h)
set(h,"color","black")
% transition zone
plot([log10(0.5), log10(300)],[0.4,0.4], "k--")
plot([log10(0.5), log10(300)],[0.35,0.35], "k--")

%% no data region
n_xxx = [log10(1.6),log10(1.6),log10(1.2), log10(1.0), log10(0.8), log10(0.8), log10(0.8),log10(0.6),log10(0.6), log10(0.5), log10(0.5), log10(0.5),log10(0.4), log10(0.4)];
n_yyy = [0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.2];
fill(n_xxx, n_yyy, [0.8,0.8,0.8],"LineStyle","-")
text(log10(1),0.3, ["   Not"; "explored"],'Rotation',0,'FontSize',12)
%%
xticks([log10(L*1000)])
xticklabels([0.4,0.5,0.6,0.8,1,1.2,1.6,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125,160,200,250,300]./5)
yticks([0.2:0.05:0.9])
yticklabels([0.2:0.05:0.9])
set(gca,'XDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
%set(gca,'YDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
xlabel("RD_{RS}")
ylabel('a/b')
box on
%% label
% text(log10(200),0.85,"EEP>0.4",'Rotation',0)
% scatter(log10(210),0.85,60,'o','MarkerEdgeColor',[1,0,0], 'MarkerFaceColor',[1,0,0])
% text(log10(200),0.8,"0.4>=EEP>0.2",'Rotation',0)
% scatter(log10(210), 0.8,60,'o','MarkerEdgeColor',[1,0.5,0], 'MarkerFaceColor',[1,0.5,0] )
% text(log10(200),0.75,"0.2>=EEP>-0.2",'Rotation',0)
% scatter(log10(210), 0.75,60,'o','MarkerEdgeColor',[0,1,0],  'MarkerFaceColor',[0,1,0])
% text(log10(200),0.7,"EEP<=-0.2",'Rotation',0)
% scatter(log10(210), 0.7,60,'o','MarkerEdgeColor',[0,0,1],  'MarkerFaceColor',[0,0,1])
%% rupture style for a/b>=0.4
text(log10(240),0.5,  ["Aseismic slip"],'Rotation',0,'FontSize',12)
% text(log10(10),0.55,"Regular seismic events",'Rotation',0)
% text(log10(63),0.45,"Symmetric-bilateral",'Rotation',40)
% text(log10(18),0.45,["Unsymmetric-";"bilateral";"and unilateral"],'Rotation',40)
% text(log10(10),0.45,"Full and partial",'Rotation',40)
% text(log10(5),0.45,"Crack-like with aftershocks ",'Rotation',40)
% text(log10(2.5),0.45,["Pulse-like with";  "aftershocks"],'Rotation',40)
% transition zone

% text(log10(30),0.375,["Transition zone of two nucleation equations"],'Rotation',0,'FontSize',12)
% rupture style for a/b<=0.35
% text(log10(200),0.25,"SSE(<0.1m/s) and Creep(<1e-8m/s)",'Rotation',40)
% text(log10(200),0.23,"Symmetric-bilateral",'Rotation',30)
% text(log10(40),0.23,["Unsymmetric-";"bilateral";"and unilateral"],'Rotation',45)
% text(log10(16),0.23,"Full and partial",'Rotation',45)
% text(log10(8),0.23,["Crack-like with"; "aftershocks"],'Rotation', 45)
% text(log10(4),0.23,["Pulse-like with";  "aftershocks"],'Rotation',45)

%% all test cases
% scatter(P_1(:,1),P_1(:,2) ,'*','MarkerEdgeColor',"#77AC30", 'MarkerFaceColor',"#77AC30")    % resolution limit
% scatter(P_1(:,1),P_1(:,2) ,'*','k')    % resolution limit
% scatter(log10(40), 0.6 ,'*','k')    % resolution limit
% scatter(log10(30), 0.65 ,'*','k')    % resolution limit
% scatter(log10(20), 0.7 ,'*','k')    % resolution limit
% scatter(log10(8), 0.8 ,'*','k')    % resolution limit
% scatter(log10(4), 0.85 ,'*','k')    % resolution limit
% scatter(log10(2), 0.9 ,'*','k')    % resolution limit
% scatter(log10(0.5), 0.9,'*','r')    % resolution limit
%% example
% scatter(log10(6), 0.8, 150, 'o','MarkerEdgeColor',[0,0,0], 'MarkerFaceColor',[1,0,0])    %  193
% scatter(log10(10), 0.7, 150,'o','MarkerEdgeColor',[0,0,0], 'MarkerFaceColor',[1,0.5,0] )    %  171
% scatter(log10(20), 0.55, 150,'o','MarkerEdgeColor',[0,0,0],  'MarkerFaceColor',[0,1,0])  %  130
% scatter(log10(20), 0.5, 150, 'o','MarkerEdgeColor',[0,0,0],  'MarkerFaceColor',[0,0,1])   %  113
lw = 1.0;
% transition
scatter(log10(6), 0.8, 180, 's','MarkerEdgeColor',[0,0,0], 'LineWidth',lw)    %  193
scatter(log10(10), 0.7, 180,'s','MarkerEdgeColor',[0,0,0],'LineWidth',lw)      %  171
scatter(log10(20), 0.55, 180,'s','MarkerEdgeColor',[0,0,0],'LineWidth',lw)  %  130
scatter(log10(20), 0.5, 180, 's','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  113
% aseismic transient
% scatter(log10(8), 0.55, 180, 'v','MarkerEdgeColor',[0,0,0])   %  167
% scatter(log10(5), 0.65, 180, 'v','MarkerEdgeColor',[0,0,0])   %  154
scatter(log10(2.5), 0.75, 200, 'p','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  128
scatter(log10(1.6), 0.8, 200, 'p','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  
scatter(log10(0.8), 0.85,200, 'p','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  
% fixed length
scatter(log10(10), 0.3, 180, 'o','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  154
scatter(log10(30), 0.3, 180, 'o','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  167
scatter(log10(80), 0.3, 180, 'o','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
% twin
scatter(log10(50), 0.35, 180, 'd','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  154
scatter(log10(63), 0.3, 180, 'd','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  167
scatter(log10(80), 0.25, 180, 'd','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
% foreshock
scatter(log10(6), 0.75, 180, '^','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  154
scatter(log10(4), 0.8, 180, '^','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  167
scatter(log10(2), 0.85, 180, '^','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
scatter(log10(1), 0.9, 180, '^','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
% fixed length with a/b>0.5
scatter(log10(8), 0.5, 180, 'h','MarkerEdgeColor',[0,0,0],'LineWidth',lw)  %  176
scatter(log10(5), 0.55, 180, 'h','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
scatter(log10(3), 0.6, 180, 'h','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
scatter(log10(2), 0.65, 180, 'h','MarkerEdgeColor',[0,0,0],'LineWidth',lw)   %  176
%%
set(gca,'TickDir', 'out')
%%
exportgraphics(gcf,'separate_nucleation_style_EEP_new.png','Resolution',600)