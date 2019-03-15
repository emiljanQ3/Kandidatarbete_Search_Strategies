function [pos_res, indice] = splitPositionData(pos_a)
    figure(3004)
    
    X = pos_a(1, 1, :);
    Y = pos_a(1, 2, :);
    start = [pos_a(1,1,1), pos_a(1,2,1)];
    slut = [pos_a(1,1,size(pos_a,3)), pos_a(1,2,size(pos_a,3))];
    
    %indice= zeros(floor(size(pos_a, 3)/10), 2); %should be enough

    X = X(:,:)';
    Y = Y(:,:)';
    plot(X,Y);
    hold on
    plot(X,Y, '.');

    plot(start(1),start(2), 'o', 'MarkerEdgeColor','r')
    plot(slut(1), slut(2),'o', 'MarkerEdgeColor','g')
    
    button = 0;
    agent = 0;
    pos_res = zeros(1,2,length(pos_a(1,1,:)))
    while(button ~= 3)
        agent = agent +1;
        [x,y,button] = ginput;
        
        scatter(x,y,'b.')
        button
        if(button == 3)
            break
        end

        dist = [sqrt((squeeze(pos_a(1,1,:)) - x(1)).^2 + squeeze((pos_a(1,2,:)) - y(1)).^2), sqrt((squeeze(pos_a(1,1,:)) - x(2)).^2 + squeeze((pos_a(1,2,:)) - y(2)).^2)];

        [~,I] = min(dist);
        
        

        
        scatter(pos_a(1, 1, I),pos_a(1, 2, I), 'r.')
        
        X = pos_a(1, 1, I(1): I(2)) ;
        Y = pos_a(1, 2, I(1):I(2));
      
        X = X(:,:)';
        Y = Y(:,:)';
        plot(X,Y, 'r');
        
        if I(2)>I(1)
            pos_res(agent,:,1:(I(2)-I(1)+1)) = pos_a(1,:,I(1):I(2));
            indice(agent, :) = [I(1); I(2)]; % start index; slut index
        else
            pos_res(agent,:,1:(I(1)-I(2)+1)) = pos_a(1,:,I(2):I(1));
            indice(agent, :) = [I(2); I(1)]; % start index; slut index
        end

    end
    indice(agent,:) = NaN
    p = 5;
end