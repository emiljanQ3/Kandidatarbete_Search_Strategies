function [pos_res, indice] = splitPositionData(pos_a, myCircle)
figure

%press start point and end point then enter to cut
%cut will turn red
%if good cut: press enter again and it turns green and is saved in indice
%if not good cut: press some other button and it turns yellow and is not
%included in indice

X = pos_a(1, 1, :);
Y = pos_a(1, 2, :);
start = [pos_a(1,1,1), pos_a(1,2,1)];
slut = [pos_a(1,1,size(pos_a,3)), pos_a(1,2,size(pos_a,3))];

X = X(:,:)';
Y = Y(:,:)';
plot(X,Y);
hold on
plot(X,Y, '.');
x = myCircle(1,1,:);
y = myCircle(1,2,:);
x = x(:,:)';
y = y(:,:)';
plot(x,y,'c')
hold on
axis equal

plot(start(1),start(2), 'o', 'MarkerEdgeColor','r')
plot(slut(1), slut(2),'o', 'MarkerEdgeColor','g')

button = 0;
agent = 0;
pos_res = zeros(1,2,length(pos_a(1,1,:)));
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
    
    if (I(1) < I(2))
        X = pos_a(1, 1, I(1): I(2)) ;
        Y = pos_a(1, 2, I(1):I(2));
    else
        X = pos_a(1, 1, I(2): I(1)) ;
        Y = pos_a(1, 2, I(2):I(1));
    end
    X = X(:,:)';
    Y = Y(:,:)';
    
    plot(X,Y, 'r');
    
    [~,~,button2] = ginput;
    
    
    if(isempty(button2)) %the space bar
        
        if I(2)>I(1)
            pos_res(agent,:,1:(I(2)-I(1)+1)) = pos_a(1,:,I(1):I(2));
            indice(agent, :) = [I(1); I(2)]; % start index; slut index
        else
            pos_res(agent,:,1:(I(1)-I(2)+1)) = pos_a(1,:,I(2):I(1));
            indice(agent, :) = [I(2); I(1)]; % start index; slut index
        end
        
        plot(X,Y, 'g');
        
    else
        plot(X,Y,'y')
        agent=agent-1;
    end
    
    
end

indice(agent,:) = [NaN,NaN];
end