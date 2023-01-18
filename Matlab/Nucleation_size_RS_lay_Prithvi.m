%% nucleation size versus half width of damage zone H
clear
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 800 300]);%×óÏÂ½ÇÎ»ÖÃ£¬¿í¸ß
%R_nuc = fzero('funx',3000)
%% r: ratio of host rock£¨shear modulus£©
%   H is the half-thickness of the low rigidity layer: from Prithvi's
%   calculation
    subplot(1,2,1)
    gamma = pi/4;  % empirical constant parameter 
    r = 0.36;   % the shear wave reduction=60%  
    mu = 3.203812032e10;  % Pa 
    mu_D = r*mu;  % Pa
    L = 8e-3;  % m
    sigma = 50e6;  % Pa
    a = 0.015;
    b = 0.019;
    h_hom = 2/pi*mu_D*L*b/sigma/(b-a)^2;   % km
    f = @(H,h_lay)(h_lay.*tanh(2*gamma.*H./h_lay+atanh(mu_D/mu)) - h_hom);
    fimplicit(f,[0 1000 2000 3800])
    xlabel('Damaged zone half-width(m)')
    ylabel('Nucleation Size(m)')
    hold on
%% 
NS = zeros(1,4); H = [75,150,300,450];
for i =1:4
    syms y
    exp = y*tanh(2*gamma*H(i)/y+atanh(mu_D/mu))-h_hom;
    y = vpasolve(exp);
    NS(i) = y;
end
scatter(H,NS,'*')
grid on
%% rr: ratio of host rock(shear wave velocity)
subplot(1,2,2)
    gamma = pi/4;  % empirical constant parameter 
    H = 150;    % half-width of fault damage zone: meters
    mu = 3.203812032e10;  % Pa 
    L = 8e-3;  % m
    sigma = 50e6;  % Pa
    a = 0.015;
    b = 0.019;
    f = @(rr,h_lay)(h_lay.*tanh(2*gamma.*H./h_lay+atanh(rr^2)) - 2/pi*rr^2*mu*L*b/sigma/(b-a)^2);
    fimplicit(f,[0.5 1 2900 3900])
    xlabel('Shear wave velocity (% of host rock)')
    ylabel('Nucleation Size(m)')
    hold on
    %% 
NS = zeros(1,4); rr = [0.6,0.7,0.8,0.9];
for i =1:4
    syms y
    exp = y*tanh(2*gamma*H/y+atanh(rr(i)^2)) - 2/pi*rr(i)^2*mu*L*b/sigma/(b-a)^2;
    y = vpasolve(exp);
    NS(i) = y;
end
scatter(rr,NS,'*')
grid on
%%
export_fig -dpng -r600 Nucleation_size_RS_lay_Prithvi


    
    