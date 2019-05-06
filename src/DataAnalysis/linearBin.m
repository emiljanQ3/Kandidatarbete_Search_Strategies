function [meanValue, binKir] = linearBin(kir,N_k, value, range)
% Makes mean in N_k different kirality bins and returns

    bins = zeros([2,N_k]);
    binWidth = (range(2)-range(1))/N_k;
    
    for i = 1:length(kir)
        binIndex = floor( (kir(i)-range(1)) / binWidth ) + 1;
        
        if binIndex > N_k || binIndex < 1
           continue 
        end
        
        bins(:, binIndex) = bins(:, binIndex) + [value(i); 1];
    end
    
    meanValue = bins(1,:)./bins(2,:);
    
    binKir = ((1:N_k) - 0.5) * binWidth + range(1);    
end