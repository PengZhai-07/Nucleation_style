% Galis,2015
% parameters
sigma_n = 1e8;   % Pa
tau_0i = 63.4e6;    % Pa
tau_0 = 56e6;   % Pa
tau_d = 0.5*sigma_n;tau_s = 0.63*sigma_n;
mu = 3e10;   % Pa   
D_c = 0.8;   % m  
%%
x = 1:0.01:5;
f = sqrt(x).*(1+(tau_0i-tau_0)/(tau_0-tau_d)*(1-sqrt(1-1./x.^2)));
plot(x,f)
f_min = min(f)
R_nuc = pi/4/f_min^2*(tau_s-tau_d)/(tau_0-tau_d)^2*mu*D_c