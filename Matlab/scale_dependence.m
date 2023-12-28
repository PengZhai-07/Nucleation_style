%%
clear
clc
close all
figure(1);
set(0,'defaultfigurecolor','w')
set(gcf,'Position',[20 20 300 250]);%左下角位置，宽高

%%
plot([log10(1) log10(1000)],[1,1],'k:', LineWidth=1)
text(log10(10), 1.2, "Lab faults")  
hold on
plot([log10(0.5) log10(500)],[2,2],'k--', LineWidth=1)
text(log10(5), 2.2, "Natural faults")  
plot([log10(0.08)  log10(60)],[3,3],'k', LineWidth=1)
text(log10(1), 3.2, "This study")  
L = [0.1,1,10,100,1000];
xticks([log10(L)])
xticklabels([0.1 1 10 100 1000])
xlabel("RD_{RS}")
box off
axis([log10(0.05) log10(1000) 0 4])
h = gca;
h.YAxis.Visible = "off";

exportgraphics(gcf,'scale_dependence.png','Resolution',600)
