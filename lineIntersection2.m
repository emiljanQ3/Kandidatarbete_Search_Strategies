function P = lineIntersection2(line1, line2)

    if (Intersects(line1(1,:),line1(2,:),line2(1,:),line2(2,:))) 
        x1 = [line1(1,1) line1(2,1)];
        y1 = [line1(1,2) line1(2,2)];
        x2 = [line2(1,1) line2(2,1)];
        y2 = [line2(1,2) line2(2,2)];
        
        [X,Y] = polyxpoly(x1',y1',x2',y2');
        P = [X,Y];
        
        return;
    end
    P = inf*0;
end