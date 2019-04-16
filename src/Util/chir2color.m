function c = chir2color(w)
    palette = getRGBGradient(0.4,0.6,1000);
    c = zeros(length(w),3);
    for i = 1:length(w)
        c(i,:) = palette(ceil((log(w(i))+1)*500),:);
    end
    
end