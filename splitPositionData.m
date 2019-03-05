function pos_res = splitPositionData(pos_a)
    figure(3004)
    
    X = pos_a(1, 1, :);
    Y = pos_a(1, 2, :);

    X = X(:,:)';
    Y = Y(:,:)';
    plot(X,Y);
    hold on
    
    button = 0;
    agent = 0;
    
    pos_res = zeros(1,2,length(pos_a(1,1,:)))
    while(button ~= 3)
        agent = agent +1;
        [x,y,button] = ginput;
        
        scatter(x,y,'b.')
        
        if(button == 3)
            break
        end

        dist = [sqrt((squeeze(pos_a(1,1,:)) - x(1)).^2 + squeeze((pos_a(1,2,:)) - y(1)).^2), sqrt((squeeze(pos_a(1,1,:)) - x(2)).^2 + squeeze((pos_a(1,2,:)) - y(2)).^2)];

        [~,I] = min(dist);

        
        scatter(pos_a(1, 1, I),pos_a(1, 2, I), 'r.')
        
        if I(2)>I(1)
            pos_res(agent,:,1:(I(2)-I(1)+1)) = pos_a(1,:,I(1):I(2));
        else
            pos_res(agent,:,1:(I(1)-I(2)+1)) = pos_a(1,:,I(2):I(1));
        end

    end
    
    p = 5;
end