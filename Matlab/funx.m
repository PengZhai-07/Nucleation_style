function fx = funx(h_lay,)
%   % kaneko 2011
% nucleation size for low-velocity layer adjacent to the fault plane and
% that of the undamaged host rock
    gamma = pi/4;  % empirical constant parameter 
    H = 500;  % thickness of low-velocity zone:m
    r = 0.8*0.8;     % 20% velocity reduction
    mu = 3.2e10;  % Pa 
    mu_D = r*mu;  % Pa
    L = 8e-3;  % m
    sigma = 50e6;  % Pa
    a = 0.015;
    b = 0.019;
    h_hom = 2/pi*mu_D*L*b/sigma/(b-a)^2   % km
    fx = h_lay*tanh(2*H*gamma/h_lay+atanh(mu_D/mu)) - h_hom;
end


