%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 600 800]);%左下角位置，宽高
%% damage zone
x = [0,1.5,1.5,0];
y = [-5,-5,15,15];
fill(x,y,"red")
hold on
%%
plot([0 0],[-5,15], 'k', "lineWidth", 2)
plot([16 16],[-5,15], 'k:', "lineWidth", 2)
plot([0 16],[15,15], 'k--', "lineWidth", 2)
plot([0 16],[-5,-5], 'k:', "lineWidth", 2)
%%
plot(2.5,0,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(2.5,0,'k.', 'MarkerSize', 20,"lineWidth", 2)
plot(-2.5,0,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(-2.5,0,'kx', 'MarkerSize', 20,"lineWidth", 2)
%%
plot([0 0],[1,13], 'r', "lineWidth", 2)
plot([-0.5 0],[13,13], 'k', "lineWidth", 2)
plot([-0.5 0],[3,3], 'k', "lineWidth", 2)
plot([0 0],[-5,15], 'k', "lineWidth", 2)
%%
text(2,5,["   Antiplane"; "Simple Shear"],"Fontsize",25)
text(-2,4,["Seismogenic"; "      Zone"],"Fontsize",25,'rotation',90)
% text(2,5,"10km","Fontsize",20,'rotation',90)
text(14,5,"D","Fontsize",25)
text(8,-6,"L","Fontsize",25)
axis equal
axis off
%%
export_fig png -r600 antiplane_model_SSA