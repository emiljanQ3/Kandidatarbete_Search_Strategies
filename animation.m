function p = animation(pos_a,obstacle,dT,colision)

    figure(2002)
    hold on
    
    L = 1;
    maxPos = max(max(pos_a,[],1),[],3);
    minPos = min(min(pos_a,[],1),[],3);
    maxmax = max(abs([maxPos minPos]));
   	plotSize = ceil((maxmax) / L);

    for i = -plotSize:plotSize
        for j = -plotSize:plotSize
            for k = 1:size(obstacle, 3)
                plot(obstacle(:,1,k)+j*L,obstacle(:,2,k)+i*L, 'k', 'LineWidth', 1)

            end
        end
    end

    for time = 2:length(pos_a(1,1,:))
        if(colision(3,time) == 1)
            scatter(colision(1,time),colision(2,time),'r')
        end
        pause(dT)
        h =  plot([pos_a(:,1,time),pos_a(:,1,time-1)]',[pos_a(:,2,time),pos_a(:,2,time-1)]','b');
        set(h, {'color'}, num2cell(jet(length(pos_a(:,1,1))),2));
    end
    
    p = 5;
end