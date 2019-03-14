function r = doublePoint(pos_a,stepSizeThreshold)
% some problems to fix
% it can get 2 nan in a row 
% handel the first and last point
% return a time vektor with the correkt time indexes if it is needed
% om vi har hål på beggä sidor om dubbelpunkten
    
    tol = 4;
    segment = pos_a;
    
    for i = 2:length(pos_a)-2
       u = pos_a(1,:,i+1)-pos_a(1,:,i);
       if(norm(u) < stepSizeThreshold)
           v2 = pos_a(1,:,i+2)-pos_a(1,:,i+1);
           v1 = pos_a(1,:,i)-pos_a(1,:,i-1);
           if(norm(v1) <tol  && norm(v2)> tol)
               segment(1,:,i+1) = inf*0; 
           elseif(norm(v1) >tol && norm(v2)< tol)
               segment(1,:,i) = inf*0; 
           else
               segment(1,:,i) = 0;
           end
       end
    end
    
    k = 0;
    for i = 1:length(segment)
        if (norm(segment(1,:,i))~=0)
            k = k+1;
            r(1,:,k) = segment(1,:,i);
        end
    end

end