%%
clear
clc
close all

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

L = [0.5,0.6,0.8,1,1.3,1.6,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125,160,200,250,300]*10^-3;

H = 0;    % half-width
NS = zeros(length(b),length(L),length(H));
Cohesive = zeros(length(b),length(L),length(H));
res_constraint = zeros(length(b),length(L),length(H));
% 
diff_NS = zeros(length(b),length(L),length(H));
EEP = zeros(length(b),length(L),length(H));

W = 5000;    % unit:m
m = 1;

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
            % get the possible values from simulations
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
                    continue
                end
            end

            fileID  = fopen(filename,'r');
            % extract EEP
            A = textscan(fileID, '%f', 4, 'HeaderLines',3);
            temp_EEP = A{1}(1);
            EEP(i,j,k) = temp_EEP;
            fclose(fileID);

        end
    end
end

%%
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 1400 800]);%左下角位置，宽高

[Y,X] = meshgrid(a_b, log10(L*1000));
A = pcolor(X,Y,EEP');
shading interp
colormap(jet)
clim([min(min(EEP)),max(max(EEP))])
% clim([-0.3, 0.3])
c = colorbar;
ylabel(c, 'EEP')


hold on 

v = [1, 2, 4, 8, 16, 28];
[c,h]=contour(X,Y,Ru',v);
clabel(c,h)
set(h,"color","black")
xticks([log10(L*1000)])
xticklabels([0.5,0.6,0.8,1,1.3,1.6,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125,160,200,250,300])
yticks([0.2:0.05:0.35,0.3781,0.4:0.05:0.9])
yticklabels([0.2:0.05:0.35,0.3781,0.4:0.05:0.9])
set(gca,'XDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
%set(gca,'YDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
xlabel('D_{c}(mm)')
ylabel('a/b')
box on
%% rupture style
text(log10(200),0.4,"SSE(<0.1m/s) and Creep(<1e-8m/s)",'Rotation',40)
text(log10(80),0.4,"Symmetric-bilateral",'Rotation',40)
text(log10(25),0.4,["Unsymmetric-";"bilateral";"and unilateral"],'Rotation',40)
text(log10(12),0.4,"Full and partial",'Rotation',40)
text(log10(6),0.4,"Crack-like with aftershocks ",'Rotation',40)
text(log10(3),0.4,["Pulse-like with";  "aftershocks"],'Rotation',40)
plot([log10(2), log10(300)],[0.3781,0.3781], "k--")
% plot([log10(0.5), log10(1000)],[0.57,0.57], "k--")

export_fig -dpng -r300 separate_nucleation_style_32_EEP
