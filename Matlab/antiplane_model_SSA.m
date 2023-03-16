%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 1000 600]);%左下角位置，宽高
%% damage zone
x = [0,2,2,0];
y = [-5,-5,15,15];
fill(x,y,"gray")
hold on
%%
plot([0 0],[-15,15], 'k', "lineWidth", 2)
plot([24 24],[-15,15], 'k:', "lineWidth", 2)
plot([0 24],[15,15], 'k--', "lineWidth", 2)
plot([0 24],[-15,-15], 'k:', "lineWidth", 2)
%%
plot(7,0,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(7,0,'k.', 'MarkerSize', 20,"lineWidth", 2)
plot(-7,0,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(-7,0,'kx', 'MarkerSize', 20,"lineWidth", 2)
%%
plot([0 0],[1,13], 'r', "lineWidth", 2)
plot([-0.5 0],[1,1], 'k', "lineWidth", 2)
plot([-0.5 0],[13,13], 'k', "lineWidth", 2)
plot([0 0],[-15,15], 'k', "lineWidth", 2)
%%
text(5,5,["   Antiplane"; "Simple Shear"],"Fontsize",20)
text(-3,2,["Seismogenic"; "      Zone"],"Fontsize",20,'rotation',90)
% text(2,5,"10km","Fontsize",20,'rotation',90)
text(19,0,"30km","Fontsize",20)
text(10,-13,"24km","Fontsize",20)
axis equal
axis off
%%
export_fig png -r600 antiplane_model_SSA