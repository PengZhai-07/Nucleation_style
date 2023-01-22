clear
clc
sigma = 50;  % MPa
f0 = 0.6;
a = 0.015;
V0 = 1e-6;   % m
b = 0.019;
Dc = 0.008;  % m
deltat = 0.1;
t = 0:0.1:4e3;
V = zeros(1,length(t));
theta = zeros(1,length(t));
tau = zeros(1,length(t));

theta(1) = 1e5;  % s
V1 = 1e-6;   % m/s
V2 = 1e-4;
V = V1;
% V*deltat/Dc    % >=1e-6 with full terms
for i
if  t(i)<200
    V = V1;
elseif 200<=t(i)<=2000
    V = V2;
else
    V = V1;
end
% before increase of V
for i = 1:length(t)
    theta(i+1) = theta(i)*exp(-V*deltat/Dc) + Dc/V*(1-exp(-V*deltat/Dc));
end
plot(t,theta(1:end-1))
for i =1:length(t)

    tau(i) =  sigma*(f0+a*log(V/V0)+b*log(V0*theta(i)/Dc));
end
plot(t,tau)



