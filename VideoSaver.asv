function VideoSaver(x,time,l1,l2)
myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)


for i=1:1:time
    j1x = l1*sin(x(1,i));
    j1y = -l1*cos(x(1,i));
    j2x = j1x + l2*sin(x(2,i)+x(1,i));
    j2y = j1y - l2*cos(x(2,i)+x(1,i));
    hold on;
    base = plot(0,0,'o','MarkerFaceColor','red','Color','blue');
    joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
    joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
    tip = plot(j2x,j2y,'*','MarkerFaceColor','red','Color','red');
    hold off;
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    delete(joint1);
    delete(joint2);
    delete(tip);
    delete(base);
    pause(0.005) %Pause and grab frame
end
close(myVideo)
end

