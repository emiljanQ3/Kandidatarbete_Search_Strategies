function [w] = getChiralitySpiral(pos_a,dT, stepSizeThreshold,n)

    str = squeeze(pos_a(1,1,:)) ~= 0;
    index = strfind(str', [1 0]);
    pos = pos_a(:,:,1:index);

    pos = doublePoint(pos_a,1);

    new_meanPath = zeros(1,2,n+1);
    count = zeros(1,n+1);
    for i = 1:size(pos, 3)-n
        i;
        %find segment
        segment = pos(:,:,i:(i+n));
        %translate to origin
        segment = segment - segment(:,:,1);
        
        %Find rotation angle
        dX = pos(:,1,i+1)-pos(:,1,i);
        dY = pos(:,2,i+1)-pos(:,2,i);
        angle = -atan2(dY,dX);

        %Rotate segment
        for t = 1:size(segment,3)
            segment(:,:,t) = [cos(angle), -sin(angle); sin(angle), cos(angle)] * segment(:,:,t)';
        end

        % idee to solve double point
        if(sqrt(dX^2+dY^2)>=1)
            new = segment;
            for k = 1:n+1
                if (~isnan(new(1,:,k)))
                    count(k) = count(k)+1;
                    new_meanPath(1,:,k) = new_meanPath(1,:,k) + new(1,:,k);
                    new_meanPath(2,:,k) = count(k);
                end
            end 
        end
    
        %Debug below
        hold on
%         if(sqrt(dX^2+dY^2)>=1)
%             figure(4)
%             hold on
         %    plot(squeeze(segment(1,1,:)),squeeze(segment(1,2,:)))
%             hold off
%             figure(3)
%             hold on
%             plot(squeeze(new(1,1,:)),squeeze(new(1,2,:)))
%             hold off
%         end

    end
    
    path(1,:) = squeeze(new_meanPath(1,1,:))./count';
    path(2,:) = squeeze(new_meanPath(1,2,:))./count';
    
    

    figure(2003)
    plot(squeeze(path(1,:)),squeeze(path(2,:)))
    axis equal
    
    T = 1:length(path);
    T = dT*(T-1);

    W = @(x) 0;
    for j = 1:length(T)
        f_x = @(x) x(3)/(x(1)^2+x(2)^2)*(x(1)-exp(-x(1)*T(j))*(x(1)*cos(x(2)*T(j))-x(2)*sin(x(2)*T(j))));
        f_y = @(x) x(3)/(x(1)^2+x(2)^2)*(x(2)-exp(-x(1)*T(j))*(x(2)*cos(x(2)*T(j))+x(1)*sin(x(2)*T(j))));
        
        w = @(x) (path(1,j)-f_x(x)).^2 + (path(2,j)-f_y(x))^2;
        W = @(x) W(x) + w(x);
        W([0.01 1 6])
    end
    
    x_0 = [0.05,1,20];
    options = optimset('MaxFunEvals',1000000,'MaxIter',1000000);
    x = fminsearch(W,x_0,options)
    
    g_x = x(3)/(x(1)^2+x(2)^2)*(x(1)-exp(-x(1)*T).*(x(1)*cos(x(2)*T)-x(2)*sin(x(2)*T)));
    g_y = x(3)/(x(1)^2+x(2)^2)*(x(2)-exp(-x(1)*T).*(x(2)*cos(x(2)*T)+x(1)*sin(x(2)*T)));

    hold on
    plot(g_x,g_y)
    
    w = x(2);
end






