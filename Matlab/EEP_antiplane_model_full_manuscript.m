%%
clear
close all
%%
figure(1) 
x = rgb2hex([0.8,0.8,0.8]);
set(0,'defaultfigurecolor',  x)
set(gcf,'Position',[20 20 1000 600]);%左下角位置，宽高

width = 16;
%%
plot([0 0],[-20,20], 'k', "lineWidth", 2)
hold on
plot([width width],[-20,20], 'k:', "lineWidth", 2)
plot([0 width],[20,20], 'k:', "lineWidth", 2)
plot([0 width],[-20,-20], 'k:', "lineWidth", 2)
%%
plot([-width -width],[-20,20], 'k:', "lineWidth", 2)
plot([0 -width],[20,20], 'k:', "lineWidth", 2)
plot([0 -width],[-20,-20], 'k:', "lineWidth", 2)
%%
plot(-3,5,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(-3,5,'k.', 'MarkerSize', 20,"lineWidth", 2)
plot(3,5,'ko', 'MarkerSize', 20,"lineWidth", 2)
plot(3,5,'kx', 'MarkerSize', 20,"lineWidth", 2)
%%
plot([0 0],[-10,10], 'r', "lineWidth", 2)
plot([-1 1],[10,10], 'k', "lineWidth", 2)
plot([-1 1],[-10,-10], 'k', "lineWidth", 2)
%%
% text(8,10,["Antiplane"; "Simple Shear"],"Fontsize",20)
text(-12,-2,["Unstable";"Asperity";"   (VW)"],"Fontsize",20,'rotation',0)
text(-5,15, "VS","FontSize", 20)
text(-5,-15, "VS", "FontSize", 20)
text(1.5,0,"W=5km","F" + ...
    "ontsize",20,'rotation',0)
% text(25,0,"10km","Fontsize",20)
% text(14,-18,"8km","Fontsize",20)
axis equal
axis off
%%
exportgraphics(gcf,'EEP_antiplane_model_manuscript.png', 'BackgroundColor',x,'Resolution',600) 