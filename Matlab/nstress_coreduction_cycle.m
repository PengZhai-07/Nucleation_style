%%
clear
close all
gamma = pi/4;  % empirical constant parameter 
r = 0.85;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
mu = 3.203812032e10;  % Pa 
H = 1500;      % m  half-width of damage zone   
%sigma = 50e6;  % Pa
a = 0.015;
b = 0.019;
L = 0.008;
multiple = [4,5,6,7];    % Pa
cos_reduction = [0.05,0.06,0.07,0.08];
W = 10000;
fid=fopen('key_par_nc.txt','w');   % rigidity ratio, half width, Lc, cohesive zone size
fprintf(fid,'ratio halfwidth L multiple cos_reduction\n');
for i = 1:length(multiple)
    for j = 1:length(cos_reduction)
         A = [r, H ,L ,multiple(i),cos_reduction(j)];
        fprintf(fid,'%.6f %.6f %.6f %.6f %.6f\n',A);
    end
end
fclose(fid);







