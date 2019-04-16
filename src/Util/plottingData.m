%% Simulated results

load('results/Final_results/2019416-1532circle_R1_t50_l0156.mat')
    T = 50;         % vid vilken tidpunkt plottar vi resultatet
    index = floor(numAreaDP*T/measurmentTime);

    figure(111)
    hold on
    c = jet(length(w));
    %c(:,2) = 0;
    for i = 1:5:size(meanAreaCovered,2)
        plot(1:numAreaDP, meanAreaCovered(:,i)./maxArea, 'color', c(i,:))
    end

    figure(112)
    hold on
    for i = 1:size(meanAreaCovered,2)
        plot(w(i), meanAreaCovered(index,i)/(maxArea),'o','color',c(i,:))
    end
    set(gca,'xscale','log')
   
    title('')


%% Experimental results
load('results/Final_results/2019416-1521_circle_medium_1agent_t45_l216.mat')
N_k = 40;


figure(110)
hold on

[~, sortOrder] = sort(abs(kir));
area_sorted = area(:,sortOrder);

color = jet(size(area_sorted,2));
for i = 1:size(area_sorted,2)
    plot(area_sorted(:,i)/maxArea,'color',color(i,:))
end

figure(109)
hold on

[meanArea,binKir] = makeMean(kir,N_k,area);

color = jet(size(meanArea,2));
for i = 1:size(meanArea,2)
    plot(meanArea(:,i)/maxArea,'color',color(i,:))
end

name = 'Experimental results';
title(name)

figure(112)
semilogx(binKir, meanArea(index,:)/maxArea,'o')
title('With mean over chirality bins against center of bin')


%%

figure(112)

axis([0.02 10 0 1.3])
title('Av agenter uppsökt area normerat mot maximal upptäckt area', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 28)

ylabel('Effektivitet', 'Interpreter', 'latex', 'fontsize', 35)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', 35)

legend('Periodisk','Begränsad','Exprimentel')



