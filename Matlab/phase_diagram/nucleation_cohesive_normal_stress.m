%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 1000 400]);%左下角位置，宽高
subplot(1,2,1)
gamma = pi/4;  % empirical constant parameter 
r = 0.8;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
mu = 3.2e10;  % Pa 
mu_D = r*mu;  % Pa
H = 500;      % m  half-width of damage zone   
L = 0.01;  % m     % change of L can affect the fault slip behavior easily!!
%sigma = 50e6;  % Pa
a = 0.015;
b = 0.019;
f = @(sigma,h_lay)(h_lay.*tanh(2*gamma.*H/h_lay+atanh(mu_D/mu))-2/pi*mu_D*L*b./sigma/(b-a)^2);   % km);
%f = @(sigma,h_lay)(h_lay.*tanh(2*gamma.*H/h_lay+atanh(mu_D/mu))-pi/4*mu_D*L/sigma/(b-a));   % km);
%f = @(sigma,h_lay)(h_lay.*tanh(2*gamma.*H/h_lay+atanh(mu_D/mu))-mu_D*L/sigma/(b-a));   % km);
fimplicit(f,[30e6 80e6 0 8000])
xlabel('Effective Normal Stress(Pa)')
ylabel('Nucleation Size(m)')
hold on
%
nn = 4;
NS = zeros(1,nn); sigma = [40e6,50e6,60e6,70e6];
for i = 1:nn
    syms y
    exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) -...
           2/pi*mu_D*L*b./sigma(i)/(b-a)^2;     % RA
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - pi*mu_D*L/4/sigma(i)/(b-a);   % RR
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - mu_D*L/sigma(i)/(b-a);   % RR
    y = vpasolve(exp);
    NS(i) = abs(y);
end
scatter(sigma,NS,'*')
W = 10000;    % width of seismogenic zone
Ru = W./NS;
grid on
%%
subplot(1,2,2)
% cohesive zone
Czone = zeros(1,nn); 
for i =1:nn
    Czone(i) = (9*pi/32)*mu*r*L/b/sigma(i);    % down limit of cohesive zone
end
sigmaa = 20e6:1e6:80e6;
Czonee = (9*pi/32)*mu*r*L/b./sigmaa;
plot(sigmaa/1e6,Czonee)
hold on
scatter(sigma/1e6,Czone,'*')
%axis([20 80 120 350])
xlabel('Effective normal stress(MPa)')
ylabel('Cohesive zone size(m)')
grid on
box on
export_fig -dpng -r600 nucleation_cohesive_normal_stress_0.019

%% calculate Ru number for different case

