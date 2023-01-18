%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 1000 800]);%左下角位置，宽高
subplot(2,1,1)
gamma = pi/4;  % empirical constant parameter 
r = 1;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
mu = 3.2e10;  % Pa 
mu_D = r*mu;  % Pa
H = 0;      % m  half-width of damage zone   
L = 0.004;  % m     % change of L can affect the fault slip behavior easily!!
%sigma = 50e6;  % Pa
a = 0.015;
sigma = 40e6;
f = @(b,h_lay)(h_lay.*tanh(2*gamma.*H/h_lay+atanh(mu_D/mu))-2/pi*mu_D*L*b./sigma/(b-a)^2);   % km);
% f = @(sigma,h_lay)(h_lay.*tanh(2*gamma.*H/h_lay+atanh(mu_D/mu))-pi/4*mu_D*L/sigma/(b-a));   % km);
% f = @(b,h_lay)(h_lay.*tanh(2*gamma.*H/h_lay+atanh(mu_D/mu))-mu_D*L/sigma/(b-a));   % km);
fimplicit(f,[0.016 0.032 0 25000])
xlabel('Evolution Effect b')
ylabel('Nucleation Size(m)')
hold on
%
b = 0.017:0.001:0.031;
nn = length(b);
NS = zeros(1,nn);
for i = 1:nn
    syms y
    exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) -...
           2/pi*mu_D*L*b(i)./sigma/(b(i)-a)^2;     % RA
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - pi*mu_D*L/4/sigma/(b(i)-a);   % RR
%   exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - mu_D*L/sigma/(b(i)-a);  %no pi
    y = vpasolve(exp);
    NS(i) = abs(y);
end
scatter(b,NS,'*')
W = 10000;    % width of seismogenic zone
Ru = W./NS;
grid on
%%
subplot(2,1,2)
% cohesive zone
Czone = zeros(1,nn); 
for i =1:nn
    Czone(i) = (9*pi/32)*mu*r*L/b(i)/sigma;    % down limit of cohesive zone
end
bb = 0.017:0.0001:0.031;
Czonee = (9*pi/32)*mu*r*L./bb/sigma;
plot(bb,Czonee)
hold on
scatter(b,Czone,'*')
%axis([20 80 120 350])
xlabel('Evolution Effect b')
ylabel('Cohesive zone size(m)')
grid on
box on
export_fig -dpng -r600 nucleation_cohesive_b_Rubin_Ampuero

%% calculate Ru number for different case

