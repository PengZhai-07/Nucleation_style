clear
close all
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 800 400]);%左下角位置，宽高
gamma = pi/4;  % empirical constant parameter 
mu = 3.203812032e10;  % Pa 
r = 0.85;
mu_D = r*mu;
sigma = 50e6;  % Pa
a = 0.015;
b = 0.019;
L = 0.008;
H = [500,1000,1500,2000,2500];
W = 12000;    % unit: m
for i = 1:5
    syms y
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) -
    %2/pi*mu_D*L*b./sigma(i)/(b-a)^2;     % RA
    %exp = y*tanh(2*gamma*H/y+atanh(mu_D/mu)) - pi*mu_D*L/4/sigma(i)/(b-a);   % RR
    exp = 1/y*tanh(2*H(i)*gamma/W*y+atanh(mu_D/mu)) -...
                       mu_D*L/sigma/(b-a)/W;    % without pi/4? 
    y = vpasolve(exp,[0,1000]);
    Ru(i) = abs(double(y));
end

