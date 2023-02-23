%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 2000 400]);%左下角位置，宽高
subplot(1,3,1)
gamma = pi/4;  % empirical constant parameter 
r = 1;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
mu = 3e10;  % Pa 
mu_D = r*mu;  % Pa
H = 0;      % m  half-width of damage zone   
sigma = 5e6;  % Pa
a= 0.005;
b = 0.01;
cell_size = 5000/2/34;   % unit: m 
asperity_criticalness = [0.1:0.1:1.0];
NS = cell_size./asperity_criticalness;
nn = length(asperity_criticalness);
Dc = zeros(1,nn);
for i = 1:nn
    syms y
    exp = NS(i)*tanh(2*gamma*H/NS(i)+atanh(mu_D/mu)) -...
           2/pi*mu_D*y*b./sigma/(b-a)^2;     % RA
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - pi*mu_D*L/4/sigma(i)/(b-a);   % RR
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - mu_D*L/sigma(i)/(b-a);   % RR
    y = vpasolve(exp);
    Dc(i) = abs(y);
end
subplot(1,3,1)
scatter(asperity_criticalness,NS,'*')
xlabel('Asperity criticalness')
ylabel('Nucleation size(m)')
hold on
grid on
subplot(1,3,2)
scatter(asperity_criticalness,Dc.*1000,'*')
xlabel('Asperity criticalness')
ylabel('Dc(mm)')
hold on
grid on
%%
subplot(1,3,3)
% cohesive zone
Czone = zeros(1,nn); 
for i =1:nn
    Czone(i) = (9*pi/32)*mu*r*Dc(i)/b/sigma;    % down limit of cohesive zone
end
scatter(asperity_criticalness,Czone,'*')
%axis([20 80 120 350])
xlabel('Asperity criticalness')
ylabel('Cohesive zone size(m)')
grid on
box on
export_fig -dpng -r600 tremor_10e6

%% calculate Ru number for different case

