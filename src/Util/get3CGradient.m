function c = get3CGradient(c1,c2,c3, trans1, trans2, l)
   

    r1 = linspace(c1(1), (c1(1)+c2(1))/2, floor(l*trans1)); 
    g1 = linspace(c1(2), (c1(2)+c2(2))/2, floor(l*trans1)); 
    b1 = linspace(c1(3), (c1(3)+c2(3))/2, floor(l*trans1));

    r2 = linspace((c1(1)+c2(1))/2, (c3(1)+c2(1))/2, floor(l*trans2)-length(r1));
    g2 = linspace((c1(2)+c2(2))/2, (c3(2)+c2(2))/2, floor(l*trans2)-length(g1)); 
    b2 = linspace((c1(3)+c2(3))/2, (c3(3)+c2(3))/2, floor(l*trans2)-length(b1));
    
    r3 = linspace((c3(1)+c2(1))/2, c3(1), l-length(r1)-length(r2));
    g3 = linspace((c3(2)+c2(2))/2, c3(2), l-length(g1)-length(g2)); 
    b3 = linspace((c3(3)+c2(3))/2, c3(3), l-length(b1)-length(b2));

    c  = [r1',g1',b1';
          r2',g2',b2';
          r3',g3',b3'];
end