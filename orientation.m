function type = orientation(p, q, r) 

    val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2)); 
  
    type = sign(val);
end 