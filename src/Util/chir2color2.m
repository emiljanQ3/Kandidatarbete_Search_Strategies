function c = chir2color2(w)
    palette = getYGBGradient(0.4,0.5,0.6,1001);
    c = zeros(length(w),3);
    for i = 1:length(w)
        if (w(i) < 0.1)
            c(i,:) = palette(1,:);
        else
            c(i,:) = palette(ceil((log10(w(i))+1)*500)+1,:);
        end
    end
    
end