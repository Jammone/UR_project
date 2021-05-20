% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear all;
close;
%% parameters
l1 = 0.5;
q1_0 = 0;
basex = 0;
time = 100;




%% plot
fig = figure();
axp = l1*cos(q1_0)+basex;
ayp = l1*sin(q1_0);

hold on
baseplot = plot(basex,0,'gs','MarkerFaceColor','red','markers',22);
axeplot = plot([basex,axp],[0,ayp],'-','MarkerFaceColor','blue','Color','blue');
xlim([-1 1]);
ylim([-1,1]);  
hold off 

pause
for k = 2:time
    figure(fig);
    delete(baseplot);
    delete(axeplot);
   
    basex = k;
    axp = l1*cos(sin(k))+basex;
    ayp = l1*sin(sin(k));

    hold on
    baseplot = plot(basex,0,'gs','MarkerFaceColor','red','markers',22);
    axeplot = plot([basex,axp],[0,ayp],'-','MarkerFaceColor','blue','Color','blue');
    xlim([-1 1]);
    ylim([-1,1]);
    hold off
    drawnow();     
end
