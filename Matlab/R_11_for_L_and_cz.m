%%
clear
close all
gamma = pi/4;  % empirical constant parameter 
r = [0.5,1/3,0.25,0.2];   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
mu = 3.2e10;  % Pa 
H = [100,200,500,1000,2000,5000];      % m  half-width of damage zone   
%sigma = 50e6;  % Pa
a = 0.015;
b = 0.019;
L = zeros(length(r),length(H));
sigma = 50e6;    % Pa
Ru = 11;
W = 10000;
fid=fopen('key_par.txt','w');   % rigidity ratio, half width, Lc, cohesive zone size
 fprintf(fid,'ratio halfwidth L co_size\n');
for i = 1:length(r)
    mu_D(i) = r(i)*mu;  % Pa
    for j = 1:length(H)
        syms y
%         exp = 1/Ru*tanh(2*H(j)*gamma/W*Ru+atanh(mu_D(i)/mu)) -...
%                2/pi*mu_D(i)*y*b/sigma/((b-a)^2)/W;     % RA
%          exp = 1/Ru*tanh(2*H(j)*gamma/W*Ru+atanh(mu_D(i)/mu)) -...
%                        pi/4*mu_D(i)*y/sigma/(b-a)/W;    % with pi/4?    % RR
          exp = 1/Ru*tanh(2*H(j)*gamma/W*Ru+atanh(mu_D(i)/mu)) -...
                       mu_D(i)*y/sigma/(b-a)/W;    % without pi/4?    % RR
         y = vpasolve(exp);
         L(i,j) = abs(y);
         CZ(i,j) = 9*pi/32 * mu_D(i)*L(i,j)/b/sigma;
         A = [r(i), H(j) ,L(i,j) , CZ(i,j)];
         fprintf(fid,'%.6f %.6f %.6f %.6f\n',A);
    end
end

fclose(fid);
