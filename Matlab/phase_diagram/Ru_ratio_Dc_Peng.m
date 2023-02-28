%% damage zone rigidity ratio and characteristic slip distance
%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 800 400]);%左下角位置，宽高
gamma = pi/4;  % empirical constant parameter 
mu = 3.203812032e10;  % Pa 
sigma = 50e6;  % Pa
a = 0.015;
b = 0.019;
r_r = linspace(5,1,21);
r = 1./r_r;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
% LL = linspace(log10(0.5),log10(125),25);  % m
% L = 10.^(LL);
L = [0.5,0.6,0.8,1,1.3,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125]*10^-3;
%L = [0.5, 1.3, 4, 12, 40, 125]*10^-3;       % m
%H = [250, 500, 1000, 1500, 2000];      % m  half-width of damage zone
H = 2000;    % half-width
NS = zeros(length(r),length(L),length(H));
W = 12000;
for i = 1:length(r)
    mu_D(i) = r(i)*mu;  % Pa
    for j = 1:length(L)   
        for k = 1:length(H)
            syms y
%             exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -...
%                2/pi*mu_D(i)*L(j)*b/sigma/(b-a)^2;
%                 exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -...
%                  pi/4*mu_D(i)*L(j)/sigma/(b-a);
               exp = 1/y*tanh(2*H(k)*gamma/W*y+atanh(mu_D(i)/mu)) -...
                       mu_D(i)*L(j)/sigma/(b-a)/W;    % without pi/4? 
             y = vpasolve(exp,[0,1000]);
            Ru(i,j,k) = y;
        end
    end
end
% Ru = W./NS;
[Y,X] = meshgrid(1./r, log10(L*1000));
% A = pcolor(X,Y,Ru');
v = [2,3,7.5,18.35,56.4,88];
%v = [2,3,5,10,15,20,30,40,50,60,70,80,100,200,400];
[c,h] = contour(X,Y,Ru',v);
xticks([log10(L*1000)])
xticklabels([0.5,0.6,0.8,1,1.3,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125])
set(gca,'XDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
set(gca,'YDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
colormap(jet)
clabel(c,h,v)
hold on 
% plot(log10(8),1.1765,'*')
xlabel('characteristic weakening distance(mm)')
ylabel('Compliancce level(G/G_{cz})')
title([num2str(H)])
export_fig -dpng -r600 Ru_ratio_Dc_Peng


