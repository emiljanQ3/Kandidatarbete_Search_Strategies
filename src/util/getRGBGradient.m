function c = getRGBGradient(trans1, trans2, l)
    c1 = [0,0,1];
    c2 = [0,1,0];
    c3 = [1,0,0];
    
    c = get3CGradient(c1,c2,c3, trans1, trans2, l);
end