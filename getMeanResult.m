%close all, clear all
k = 5
c = load('results/Lab/Tefat_c1agent.txt');
[kir, I] = sort(c(:,1));

A = c((I),2)

kir_m = movmean(kir,k);
A_m = movmean(A,k);


figure
semilogx(abs(kir),A,'o')
title('no mean')
figure
semilogx(abs(kir_m), A_m, 'o')
title('with mean')

axis([0.01 10 0 1])

      


