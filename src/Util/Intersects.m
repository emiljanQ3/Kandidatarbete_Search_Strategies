function bool = Intersects(p1, q1, p2, q2)
    o1 = orientation(p1, q1, p2); 
    o2 = orientation(p1, q1, q2); 
    o3 = orientation(p2, q2, p1); 
    o4 = orientation(p2, q2, q1);
    
    o4;
    if (o1 ~= o2 && o3 ~= o4) 
        bool = true;
        return;
    end
    bool = false;
end
