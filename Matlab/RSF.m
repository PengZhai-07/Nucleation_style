%%
clear
clc
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 600 300]);%左下角位置，宽高

%%
sigma = 50;  % MPa
f0 = 0.6;
a = 0.015;
V0 = 1e-6;   % m
Dc = 1e-3;  % m
deltat = 0.01;
t = 0:deltat:10000;     % s
t1 = 500;
t2 = 600;
V = zeros(1,length(t));
D = zeros(1,length(t));
theta = zeros(1,length(t));
tau = zeros(1,length(t));
V1 = 1e-6;   % m/s
V2 = 1e-4;
theta(1) = Dc/V1;  % initial steady state 
% V*deltat/Dc    % >=1e-6 with full terms
for i = 1:length(t)
    if  t(i)<=t1 || t(i)>=t2
        V(i) = V1;
    else
        V(i) = V2;
    end
end

%% weakening
b = 0.025;
for i = 1:length(t)
    theta(i+1) = theta(i)*exp(-V(i)*deltat/Dc) + Dc/V(i)*(1-exp(-V(i)*deltat/Dc));
end
% plot(t,theta(1:end-1))
for i =1:length(t)
    tau(i) =  sigma*(f0+a*log(V(i)/V0)+b*log(V0*theta(i)/Dc));
end
for i =1:length(t)
    if t(i)<=t1
        D(i) = t(i).*V1;
    elseif t(i)>t1&&t(i)<t2
        D(i) = t1*V1+V2*(t(i)-t1);
    else
        D(i) = t1*V1+(t2-t1)*V2+(t(i)-t2)*V1;
    end
end
plot(D,tau/sigma,'r')
hold on

%% strengthening
b = 0.010;
for i = 1:length(t)
    theta(i+1) = theta(i)*exp(-V(i)*deltat/Dc) + Dc/V(i)*(1-exp(-V(i)*deltat/Dc));
end
% plot(t,theta(1:end-1))
for i =1:length(t)
    tau(i) =  sigma*(f0+a*log(V(i)/V0)+b*log(V0*theta(i)/Dc));
end
for i =1:length(t)
    if t(i)<=t1
        D(i) = t(i).*V1;
    elseif t(i)>t1&&t(i)<t2
        D(i) = t1*V1+V2*(t(i)-t1);
    else
        D(i) = t1*V1+(t2-t1)*V2+(t(i)-t2)*V1;
    end
end
plot(D,tau/sigma,'b--')
%%
legend('Velocity weakening(a-b<0)','Velocity strengthening(a-b>0)')
xlabel('Loading distance(m)')
ylabel('Friction coefficient')
title("Velocity change 10^{-6}(m/s)-10^{-4}(m/s)-10^{-6}(m/s)")
export_fig -dpng -r600 RSF