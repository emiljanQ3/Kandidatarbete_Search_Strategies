function [meanArea, binKir, limit] = makeMean(kir,N_k, area)
% Makes mean in N_k different kirality bins and returns

    N = size(area,1);

    Mi = -1;
    Ma = max(log10(abs(kir)));
    L = (Ma-Mi)/N_k;
    limit=zeros(N_k,1);
    limit(1)=10
    meanArea = zeros(N,N_k);
    sumKir = zeros(1,N_k);
    binKir = zeros(1,N_k);
    count = zeros(1,N_k);
    for i = 1:N_k
        limit(i) = 10^(Mi + (i-1)*L);
        for j = 1:length(kir)
            if ((i-1)*L + Mi <= log10(abs(kir(j))) && log10(abs(kir(j))) <= Mi + (i)*L )
                meanArea(:,i) = meanArea(:,i) + area(:,j);
                sumKir(i) = sumKir(i) +abs(kir(j));
                count(i) = count(i) +1;
            end
        end
        meanArea(:,i) = meanArea(:,i)/count(i);
        binKir(i) = 10^(Mi + L*(i+1/2));
    end
    limit(i+1) = 10^Ma;
    
end