function resultPos = moveAgentPrime(pos,newPos,obstacle,L, threshold)
        
        if(norm(pos-newPos) < threshold)
            resultPos = pos;
            return
        end
        
        intersectionDist = inf; 
         p_i = NaN; %intersection point
         wall = zeros(2,2);
         
         for i = 1:size(obstacle,3)
             point = lineIntersection([newPos,pos], obstacle(:,:,i)*L);
             if ~isnan(point) 
                 pointDist = norm(pos-point);
                 if pointDist < intersectionDist
                     intersectionDist = pointDist;
                     p_i = point;
                     wall = obstacle(:,:,i);
                 end
             end
         end
          
         if(isnan(p_i))
             resultPos = newPos;
            return 
         end
         
         tangent = wall(:,2) - wall(:,1);
         tangent = tangent/norm(tangent);
         normal = [0,-1;1,0]*tangent;
         
         newPos = p_i + tangent*dot((pos-p_i), tangent) - 0.00001*normal*dot((pos-p_i), normal);
         
         resultPos = moveAgentPrime(pos, newPos, obstacle, L, threshold);
end
