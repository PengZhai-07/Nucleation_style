%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 1000 600]);%左下角位置，宽高
gamma = pi/4;  % empirical constant parameter 
mu = 3.203812032e10;  % Pa 
sigma = 50e6;  % Pa
a = 0.015;
b = 0.019;
r_r = linspace(5,1,21);
r = 1./r_r;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
mu_D = zeros(1,length(r));
% LL = linspace(log10(0.5),log10(125),25);  % m
% L = 10.^(LL);
% L = [0.5,0.6,0.8,1,1.3,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125]*10^-3;
%L = [0.5, 1.3, 4, 12, 40, 125]*10^-3;       % m
H = [10:2:40,50:10:90,100:20:400,500:100:900,1000:200:4000,5000:1000:1e4];  % m
L = zeros(length(r),length(H));
W = 10000; Ru = 11;
NS = W/Ru;
syms y
for i = 1:length(r)
    mu_D(i) = r(i)*mu;  % Pa
%     for j = 1:length(L)   
        for k = 1:length(H)
%             exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -...
%                2/pi*mu_D(i)*L(j)*b/sigma/(b-a)^2;
%                  exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -...
%                    pi/4*mu_D(i)*L(j)/sigma/(b-a);
               exp = 1/Ru*tanh(2*H(k)*gamma/W*Ru+atanh(mu_D(i)/mu)) -...
                       mu_D(i)*y/sigma/(b-a)/W;    % without pi/4? 
             S = vpasolve(exp,y,[0,1000]);
             L(i,k) = S;
        end
end
% Ru = W./NS;
[Y,X] = meshgrid(1./r, log10(H));
% A = pcolor(X,Y,Ru');
v = [7,8,10,12,16,20,25,30];
%v = [2,3,5,10,15,20,30,40,50,60,70,80,100,200,400];
[c,h] = contour(X,Y,1e3.*L',v);
colormap(jet)
clabel(c,h,v)
H = [10:10:90,100:100:900,1000:1000:1e4];  % m
set(gca,'XTick',log10(H));
xticklabels(["","","","","","","","","","10^{2}","","","","5*10^{2}","","","","",...
    "10^{3}","","","","5*10^{3}","","","","","10^{4}"])
xlabel('Fault zone half-width(m)')
ylabel('Compliancce level(G/G_{cz})')
hold on 
plot([log10(NS),log10(NS)],[1,5],':r','LineWidth',1)
title(['R_{u}=',num2str(Ru)])
grid on
export_fig -dpng -r600 Ru_11_Peng

