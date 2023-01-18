%%
clear
clc
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 800 400]);%左下角位置，宽高
%%
gamma = pi/4;  % empirical constant parameter 
mu = 3.20e10;  % Pa 
sigma = 40e6;  % Pa
a = 0.015;
a_b = 0.5:0.05:0.95;
b = a./a_b;
r = 1;   % the shear wave reduction=20%  1-0.2=0.8  r is the rigidity ratio
% LL = linspace(log10(0.5),log10(125),25);  % m
% L = 10.^(LL);
L = [0.5,0.6,0.8,1,1.3,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125]*10^-3;
%L = [0.5, 1.3, 4, 12, 40, 125]*10^-3;       % m
%H = [250, 500, 1000, 1500, 2000];      % m  half-width of damage zone
H = 0;    % half-width
NS = zeros(length(r),length(L),length(H));
W = 10000;
m = 0;
for i = 1:length(b)
    mu_D = mu;  % Pa
    for j = 1:length(L)   
        for k = 1:length(H)
            syms y
%             exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -...
%                2/pi*mu_D(i)*L(j)*b/sigma/(b-a)^2;
%                 exp = y*tanh(2*gamma*H(k)/y+atanh(mu_D(i)/mu)) -...
%                  pi/4*mu_D(i)*L(j)/sigma/(b-a);
               exp = 1/y*tanh(2*H(k)*gamma/W*y+atanh(mu_D/mu)) -...
                       mu_D*L(j)/sigma/(a/a_b(i)-a)/W;    % without pi/4? 
            y = double(vpasolve(exp,[0,1000]));
            if (3 <= y) && (y<= 18.35)
                m = m+1;
                P(m,:) = [log10(L(j)*1000), a_b(i), L(j),b(i)];
            end
            Ru(i,j,k) = y;
%            Ru(i,j,k) = W/(mu_D*L(j)/sigma/(b(i)-a));
        end
    end
end
% Ru = W./NS;
[Y,X] = meshgrid(a_b, log10(L*1000));
% A = pcolor(X,Y,Ru');
v = [2,3,7.5,18.35,56.4,88];
% v = [2,3,5,10,15,20,30,40,50,60,70,80,100,200,400];
[c,h] = contour(X,Y,Ru',v);
xticks([log10(L*1000)])
xticklabels([0.5,0.6,0.8,1,1.3,1.5,2,2.5,3,4,5,6,8,10,12,16,20,25,30,40,50,63,80,100,125])
set(gca,'XDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
% set(gca,'YDir','reverse');        %将x轴方向设置为反向(从右到左递增)。
colormap(jet)
clabel(c,h,v)
xlabel('Characteristic weakening distance(mm)')
ylabel('a/b')
box on
hold on

%% experment points
scatter(P(:,1),P(:,2) ,'*','r' )
%%
% yy = a./[0.018,0.019,0.020,0.021,0.022,0.023,0.025,0.027,0.029, 0.031];
% xx = log10(10)*ones(1,length(yy));
% for i = 1:length(yy)
%     if  yy(i) > a/0.022    
%         if yy(i) > a/0.021
%             scatter(xx(i),yy(i),'*','r')      % bilateral and expanding
%         else
%             scatter(xx(i),yy(i),'*','b')      % bilateral and fixed length
%         end
%     elseif yy(i) > a/0.030
%          scatter(xx(i),yy(i),'o','b')      % unilateral and fixed length
%     else
%         scatter(xx(i),yy(i),'^','b')      % partial and fixed length
%     end
% end
% xxx_1 = [3, 4, 6, 8, 10, 12, 16, 20, 25];
% yy_1 = a./0.021*ones(1,length(xxx_1));
% xx_1 = log10(xxx_1);
% for i = 1:length(xx_1)
%     if  xxx_1(i) > 6
%         if  xxx_1(i) >10 || xxx_1(i)==8
%             scatter(xx_1(i),yy_1(i),'*','r')  % bilateral and expanding
%         else
%             scatter(xx_1(i),yy_1(i),'*','b')  % bilateral and fixed length
%         end
%     elseif xxx_1(i) > 3
%         scatter(xx_1(i),yy_1(i),'o','b')  % unilateral and fixed length
%     else
%         scatter(xx_1(i),yy_1(i),'^','b')  % partial and fixed length
%     end
% end
% xxx_2 = [2, 2.5, 3, 4, 6, 8, 10, 12, 16];
% yy_2 = a./0.019*ones(1,length(xxx_2));
% xx_2 = log10(xxx_2);
% for i = 1:length(xx_2)
%     if  xxx_2(i) > 4   
%         scatter(xx_2(i),yy_2(i),'*','r')  % bilateral and expanding
%     elseif xxx_2(i) == 4    
%         scatter(xx_2(i),yy_2(i),'o','r')  % unilateral and expanding
%     elseif xxx_2(i) > 2            
%         scatter(xx_2(i),yy_2(i),'o','b')  % unilateral and fixed length
%     else
%         scatter(xx_2(i),yy_2(i),'^','b')  % partial and fixed length 
%     end
% end

% title([num2str(H)])
%% 
% export_fig -dpng -r600 Nucleation_size_phase_diagram_b_L_universal

%% output the bash script for sbatch in Great lakes

%%
fid  = fopen('./whole_space.sh','wt');%
fprintf(fid,'#!/bin/bash\n\n');    
[u, v] = size(P);
N = 133;
NN = 2; % 9 task per node: 4*9=36 
nn = 4;  % number of processors for each node

fprintf(fid,['#SBATCH --array=',num2str(N),'-',num2str(N+NN-1),'\n']);
fprintf(fid,'#SBATCH --nodes=1\n');   % only one node because current MPI doesn't work
fprintf(fid,'#SBATCH --mem=20000m\n'); 
fprintf(fid,'#SBATCH --time=14-00:00:00\n');
fprintf(fid,'#SBATCH --partition=standard\n');
fprintf(fid,'#SBATCH --account=yiheh1\n');   
fprintf(fid,['#SBATCH --ntasks-per-node=',num2str(nn),'\n']);    % multitask(openmp)
fprintf(fid,'#SBATCH --job-name=case');
fprintf(fid,'#SBATCH --output=/home/%%u/log/%%x-%%j.log\n');
fprintf(fid,'#SBATCH --error=/home/%%u/log/error-%%x-%%j.log\n\n');


fprintf(fid,['srun julia --threads ',num2str(nn),' run.jl 0.8 500 ',num2str(P(i,3)),' 4 0.00 ',num2str(P(i,4)),'\n']);
fprintf(fid,['srun julia --threads ',num2str(nn),' run.jl 0.8 500 ',num2str(P(i,3)),' 4 0.00 ',num2str(P(i,4)),'\n']);
   
