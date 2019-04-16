function r = doublePoint(pos_a,dpTol,gapTol)
% some problems to fix
% it can get 2 nan in a row 
% handel the first and last point
% om vi har hål på beggä sidor om dubbelpunkten
    

    %Remove double points
    N_zeros = 0;
    for i = 1:length(pos_a)-1
       u = pos_a(1,:,i+1)-pos_a(1,:,i);
       if(norm(u) < dpTol)
            pos_a(1,:,i) = 0;
            N_zeros = N_zeros + 1;
       end
    end
    
    %remove zeros
    r = zeros(1,2,length(pos_a)-N_zeros-1);
    k = 0;
    for i = 1:length(pos_a)
        if (norm(pos_a(1,:,i))~=0)
            k = k+1;
            r(1,:,k) = pos_a(1,:,i);
        end
    end
    
    %Interpolate in gaps
    results = zeros(1,2,2*length(r));
    k = 1;
    for j = 1:length(r)-1
        
           results(1,:,k) = r(1,:,j);
           k = k+1;
           
           u = r(1,:,j+1)-r(1,:,j);
           if(norm(u) > gapTol)
               results(1,:,k) = (r(1,:,j+1)+r(1,:,j))/2;
               k = k+1;
           end
    end
    
    %remove zeros
    k = 0;
    for i = 1:length(results)
        if (norm(results(1,:,i))~=0)
            k = k+1;
            r(1,:,k) = results(1,:,i);
        end
    end

end