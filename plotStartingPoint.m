function [] = plotStartingPoint(startInfo, obstacle, v)

    hold on
    %PlotObstacle
    for k = 1:size(obstacle, 3)
        plot(obstacle(:,1,k),obstacle(:,2,k), 'k', 'LineWidth', 1)
    end
    
    %Plot start point
    w = startInfo(1);
    rot = startInfo(2);
    x = startInfo(3);
    y = startInfo(4);
    d = v/w/pi;
    scatter(x+0.05*cos(rot),y+0.05*sin(rot), 'r.')
    scatter(x,y, 'bo')
    
    title(join(["w: ", string(w), "   d: ", string(d), "   rot: ", string(rot), "   X: ", string(x),  "  Y: ", string(y)]));
    xlim([0,1]);
    ylim([0,1]);
end