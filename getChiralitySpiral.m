function [w] = getChiralitySpiral(pos_a,dT, stepSizeThreshold,tol,n)

    pos = doublePoint(pos_a,stepSizeThreshold,tol);
    pathSum = zeros(1,2,n);
    for i = 1:size(pos, 3)-(n-1)
        %find segment
        segment = pos(:,:,i:(i+n-1));
        %translate to origin
        segment = segment - segment(:,:,1);
        
        %Find rotation angle
        u = pos(:,:,i+1)-pos(:,:,i);
        angle = -atan2(u(2),u(1));

        %Rotate segment
        for t = 1:size(segment,3)
            segment(:,:,t) = [cos(angle), -sin(angle); sin(angle), cos(angle)] * segment(:,:,t)';
        end

        pathSum = pathSum + segment(1,:,:);

        %Debug below
   
%              figure(4)
%              hold on
%             plot(squeeze(segment(1,1,:)),squeeze(segment(1,2,:)))
%              hold off



    end

    meanPath(1,:) = squeeze(pathSum(1,1,:))./i;
    meanPath(2,:) = squeeze(pathSum(1,2,:))./i;


    figure(2003)
    plot(squeeze(meanPath(1,:)),squeeze(meanPath(2,:)))
    axis equal
    
    T = 1:length(meanPath);
    T = dT*(T-1)

    W = @(x) 0;
    for j = 1:length(T)
        f_x = @(x) x(3)/(x(1)^2+x(2)^2)*(x(1)-exp(-x(1)*T(j))*(x(1)*cos(x(2)*T(j))-x(2)*sin(x(2)*T(j))));
        f_y = @(x) x(3)/(x(1)^2+x(2)^2)*(x(2)-exp(-x(1)*T(j))*(x(2)*cos(x(2)*T(j))+x(1)*sin(x(2)*T(j))));
        
        w = @(x) (meanPath(1,j)-f_x(x)).^2 + (meanPath(2,j)-f_y(x)).^2;
        W = @(x) W(x) + w(x);
    end
    
    guess = getCirality(pos,dT,stepSizeThreshold);
    x_0 = [0.05,guess,100];
    options = optimset('MaxFunEvals',1000000,'MaxIter',1000000);
    x = fminsearch(W,x_0,options)
    W(x)/(length(T)^2)
    g_x = x(3)/(x(1)^2+x(2)^2)*(x(1)-exp(-x(1)*T).*(x(1)*cos(x(2)*T)-x(2)*sin(x(2)*T)));
    g_y = x(3)/(x(1)^2+x(2)^2)*(x(2)-exp(-x(1)*T).*(x(2)*cos(x(2)*T)+x(1)*sin(x(2)*T)));

    hold on
    plot(g_x,g_y)
   
    w = x(2);
end






