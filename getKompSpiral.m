function [w,D_r,v] = getKompSpiral(pos_a,dT, stepSizeThreshold,tol,~)
    count = 0;
    for agent = 1:size(pos_a,1)
%         figure(9112)
%         hold on
%         plot(squeeze(pos_a(agent,1,:)),squeeze(pos_a(agent,2,:)),'.')
%         axis equal
        temp = doublePoint(pos_a(agent,:,:),stepSizeThreshold,tol);
        count = count + length(temp);
%         figure(9111)
%         hold on
%         plot(squeeze(temp(1,1,:)),squeeze(temp(1,2,:)),'.')
%         axis equal
    end
    n = ceil(count/100);

    if n<=3
       n = 3;
    end 
    
    count = zeros(1,n);
    count(1) = 1;
    pathSum = zeros(2,n);
    for p = 2:n
        for agent = 1:size(pos_a,1)
            pos = doublePoint(pos_a(agent,:,:),stepSizeThreshold,tol);
            for i = 1:size(pos, 3)-(p-1)
                segment = squeeze(pos(1,:,i:(i+p-1)));
                segment = segment - segment(:,1);

                %Rotat segment
                u = pos(1,:,i+1)-pos(1,:,i);
                angle = -atan2(u(2),u(1));
                if i > 1
                    v = pos(1,:,i+1)-pos(1,:,i-1);
                    angle = -atan2(v(2),v(1));
                end
                for t = 1:size(segment,2)
                    segment(:,t) = [cos(angle), -sin(angle); sin(angle), cos(angle)] * segment(:,t);
                end

                count(p) = count(p) + 1;
                pathSum(:,p) = pathSum(:,p) + segment(:,p);

%                 figure(4)
%                 hold on
%                 plot(squeeze(segment(1,:)),squeeze(segment(2,:)))
%                 hold off
            end
        end
    end
    meanPath(1,:) = pathSum(1,:)./count;
    meanPath(2,:) = pathSum(2,:)./count;
  
    T = 1:length(meanPath);
    T = dT*(T-1);

    W = @(x) 0;
    for j = 1:length(T)
        f_x = @(x) x(3)/(x(1)^2+x(2)^2)*(x(1)-exp(-x(1)*T(j))*(x(1)*cos(x(2)*T(j))-x(2)*sin(x(2)*T(j))));
        f_y = @(x) x(3)/(x(1)^2+x(2)^2)*(x(2)-exp(-x(1)*T(j))*(x(2)*cos(x(2)*T(j))+x(1)*sin(x(2)*T(j))));
        
        w = @(x) (meanPath(1,j)-f_x(x)).^2 + (meanPath(2,j)-f_y(x)).^2;
        W = @(x) W(x) + w(x)/j;
    end
    
    guess = getComplexCirality(pos_a,dT,stepSizeThreshold);
    x_0 = [0.05,guess,100];
    options = optimset('MaxFunEvals',1000000,'MaxIter',1000000);
    x = fminsearch(W,x_0,options);
   % W(x)/(length(T)^2);
    g_x = x(3)/(x(1)^2+x(2)^2)*(x(1)-exp(-x(1)*T).*(x(1)*cos(x(2)*T)-x(2)*sin(x(2)*T)));
    g_y = x(3)/(x(1)^2+x(2)^2)*(x(2)-exp(-x(1)*T).*(x(2)*cos(x(2)*T)+x(1)*sin(x(2)*T)));

    
    figure(2003)
    hold on
    plot(meanPath(1,:),meanPath(2,:))
    axis equal
    plot(g_x,g_y)
    

    w = x(2);
    D_r = x(1);
    v = x(3);
end

