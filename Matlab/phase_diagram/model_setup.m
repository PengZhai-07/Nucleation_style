%%
clear
close all
figure(1)
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 1000 600]);%左下角位置，宽高
%%
plot([0 0],[-15,15], 'k', "lineWidth", 2)
hold on
plot([24 24],[-15,15], 'k:', "lineWidth", 2)
plot([0 24],[15,15], 'k:', "lineWidth", 2)
plot([0 24],[-15,-15], 'k:', "lineWidth", 2)
%%
plot(7,0,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(7,0,'k.', 'MarkerSize', 20,"lineWidth", 2)
plot(-7,0,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(-7,0,'kx', 'MarkerSize', 20,"lineWidth", 2)
%%
plot([0 0],[-5,5], 'r', "lineWidth", 2)
plot([-1 1],[5,5], 'k', "lineWidth", 2)
plot([-1 1],[-5,-5], 'k', "lineWidth", 2)
%%
text(5,5,["Antiplane"; "Simple Shear"],"Fontsize",20)
text(-3,-3,["Unstable"; "Asperity"],"Fontsize",20,'rotation',90)
text(2,-2,"5km","Fontsize",20,'rotation',90)
text(19,0,"10km","Fontsize",20)
text(10,-13,"8km","Fontsize",20)
axis equal
axis off
%%
export_fig png -r600 model_setup