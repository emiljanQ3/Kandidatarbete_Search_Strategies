function p = animateExperiment(M,obstacle,dT)

    maxPos = max(max(M,[],1),[],3);
    minPos = min(min(M,[],1),[],3);

    figure(2002)
    hold on
    axis([minPos(1) maxPos(1) minPos(2) maxPos(2)]);
    
    % Add plot of obstacles

    for time = 2:length(M(1,1,:))
        pause(dT)
        h =  plot([M(:,1,time),M(:,1,time-1)]',[M(:,2,time),M(:,2,time-1)]','b');
        set(h, {'color'}, num2cell(jet(length(M(:,1,1))),2));
    end
    
    p = 5;
end