function [w_average, v,k] = getCirality(pos_a,dT, stepSizeThreshold)
  
    n = 0;
    t = 1;
    
    v = zeros(1,length(pos_a)-1);            %INIT
    T = zeros(1,length(pos_a)-1);            %INIT
    orientation = zeros(1,length(pos_a)-1);  %INIT
    
    
    u = pos_a(1,:,2) - pos_a(1,:,1);
    v(1) = norm(u)/dT;
    orientation(1) = atan2(u(2),u(1));
    T(1) = 0;
    
    if(norm(u)<stepSizeThreshold)
        orientation(1) = 0;
    end
    
    for i = 2:length(pos_a)-1
        
        u = pos_a(1,:,i+1) - pos_a(1,:,i);
        v(i) = norm(u)/dT;
        
        if (norm(u) > stepSizeThreshold)
            t = t+1;
            T(t) = dT*(i-1);
            angle =  atan2(u(2),u(1));
            if((orientation(t-1)-(angle+2*pi*n)) > 5)
                n = n+1;
            elseif((orientation(t-1) - (angle+2*pi*n)) < -5)
                n = n-1;
            end
            orientation(t) = angle + 2*n*pi;
        end
        
    end
    
    T = T(1:t);
    orientation = orientation(1:t);
    
    %dw = zeros(1,length(orientation)-1); %INIT 
    k = 0;
    for i = 1:length(orientation)-1
        
        if(orientation(i+1)~= 0 && orientation(i)~= 0)  
            k = k+1;
            dw(k) = (orientation(i+1)-orientation(i))/(T(i+1)-T(i));  
        end
        
    end

    w_average = mean(dw);
    v = mean(v);
    
end