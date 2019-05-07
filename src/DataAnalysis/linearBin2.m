function [meanValue, binKir1, binKir2] = linearBin2(kir1, kir2, N_k, value, range)
% Makes mean in N_k different kirality bins and returns

    bins = zeros([N_k, N_k, 2]);
    binWidth = (range(2)-range(1))/N_k;
    
    for i = 1:length(kir1)
        binIndex = floor( ([kir1(i), kir2(i)]-range(1)) / binWidth ) + 1;
        
        if max(binIndex) > N_k || min(binIndex) < 1
           continue 
        end
        
        bins(binIndex(1), binIndex(2), 1) = bins(binIndex(1), binIndex(2), 1) + value(i);
        bins(binIndex(1), binIndex(2), 2) = bins(binIndex(1), binIndex(2), 2) + 1;
    end
    
    meanValue = bins(:,:,1)./bins(:,:,2);
    
    binKir1 = ((1:N_k) - 0.5) * binWidth + range(1);
    binKir2 = binKir1;
    
    
end