%%
clear
clc
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 800 500]);%左下角位置，宽高

%%
sigma = 50;  % MPa
f0 = 0.6;
a = 0.015;
V0 = 1e-6;   % m
Dc = 1e-2;  % m
deltat = 0.001;
t = 0:deltat:80000;     % s
t1 = 5000;
t2 = 12000;
V = zeros(1,length(t));
D = zeros(1,length(t));
theta = zeros(1,length(t));
tau = zeros(1,length(t));
V1 = 1e-6;   % m/s
V2 = 1e-5;
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
b = 0.020;
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
plot(D,tau/sigma,'r',LineWidth=1 )
hold on
plot([-0.02 0.0],[0.6,0.6],'r', LineWidth=1)

%% strengthening
% b = 0.010;
% for i = 1:length(t)
%     theta(i+1) = theta(i)*exp(-V(i)*deltat/Dc) + Dc/V(i)*(1-exp(-V(i)*deltat/Dc));
% end
% % plot(t,theta(1:end-1))
% for i =1:length(t)
%     tau(i) =  sigma*(f0+a*log(V(i)/V0)+b*log(V0*theta(i)/Dc));
% end
% for i =1:length(t)
%     if t(i)<=t1
%         D(i) = t(i).*V1;
%     elseif t(i)>t1&&t(i)<t2
%         D(i) = t1*V1+V2*(t(i)-t1);
%     else
%         D(i) = t1*V1+(t2-t1)*V2+(t(i)-t2)*V1;
%     end
% end
% plot(D,tau/sigma,'b--')

%%
plot([0.005,0.015], [0.595, 0.595], 'k-')
plot([0.005,0.08], [0.6+a*log(10), 0.6+a*log(10)], 'k--')
plot([0.005,0.015], [0.595, 0.595], 'k-')
annotation("arrow", [0.095/0.18, 0.095/0.18], [0.75, (0.045+a*log(10))/0.091])
annotation("arrow", [0.095/0.18, 0.095/0.18], [0.6, (0.054+a*log(10)-b*log(10))/0.091])
text(0.007,0.59,'D_{RS}')
text(-0.01,0.62,'$$aln(\frac{V_{2}}{V_{1}})$$', Interpreter='latex')
text(0.065,0.612,'$$bln(\frac{V_{2}}{V_{1}})$$', Interpreter='latex')
text(-0.015,0.57,'SLOW(V_{1})')
text(0.03,0.57,'FAST(V_{2})')
text(0.1,0.57,'SLOW(V_{1})')
%%
% legend('Velocity weakening(a-b<0)','Velocity strengthening(a-b>0)', ['D_{RS}=', num2str(Dc),'m'])
legend('Velocity weakening(a<b)')
legend boxoff
xlabel('Loading distance(m)')
ylabel('Friction coefficient')
% title("Velocity stepping 10^{-6}(m/s)-10^{-5}(m/s)-10^{-6}(m/s)")
title("Velocity stepping experiments (V_{1}<V_{2})")
box off
axis off
exportgraphics(gcf,'RSF.png','Resolution',600)