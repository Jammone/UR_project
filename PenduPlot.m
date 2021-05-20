function  PenduPlot(x,time,l1,l2)
%PLOT Summary of this function goes here
%   Detailed explanation goes here
    fig = figure();
    j1x = l1*cos(x(1,1));
    j1y = l1*sin(x(1,1));
    j2x = j1x + l2*cos(x(2,1)+x(1,1));
    j2y = j1y + l2*sin(x(2,1)+x(1,1));
    
    hold on
    joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
    joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
    tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
    xlim([-1 1]);
    ylim([-1,1]); 
    hold off 

    for k = 1:time
        figure(fig)
        delete(joint1);
        delete(joint2);
        delete(tip);
        j1x = l1*cos(x(1,k));
        j1y = l1*sin(x(1,k));
        j2x = j1x + l2*cos(x(2,k)+x(1,k));
        j2y = j1y + l2*sin(x(2,k)+x(1,k));

        hold on
        joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
        joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
        tip = plot(j2x,j2y,'*','MarkerFaceColor','red','Color','red');
        xlim([-1 1]);
        ylim([-1,1]);
        hold off
        drawnow();
    
    end
end

