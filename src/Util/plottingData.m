%% Simulated results

load('results/Final_results/2019417-1033_circle_R1_t50_l0156.mat')

    T       = 50;         % vid vilken tidpunkt plottar vi resultatet mot kiralitet
    T2      = 50;         % Hur lång tid plottar vi arean 
    index2  = floor(numAreaDP*T2/measurmentTime);           
    index   = floor(numAreaDP*T/measurmentTime);

    figure(111)
    hold on
   % c = chir2color(abs(w(1:end)));
    c = getRGBGradient(0.4,0.5,0.6,100)
    map = colormap(c);
    for i = 1:size(meanAreaCovered,2)
        plot(1:numAreaDP, meanAreaCovered(1:index2,i)./maxArea, 'color', c(i,:))
    end
    bar = colorbar
    bar.Ruler.Scale = 'log';
    bar.Ruler.MinorTick = 'on';
    caxis([0.1 10])

    figure(112)
    hold on
%     for i = 1:size(meanAreaCovered,2)
%         plot(w(i), meanAreaCovered(index,i)/(maxArea),'0','color',c(i,:))
%     end
    %plot(w,meanAreaCovered(index,i)/maxArea,'b','Linewidth',3)
    for i = 1:size(meanAreaCovered,2)-1
        
        plot([w(i+1) w(i)],[meanAreaCovered(index,i+1), meanAreaCovered(index,i)]/maxArea,'color',c(i,:),'Linewidth',3)
    end
    set(gca,'xscale','log')
    axis([0.1 10 0 1])
    title('')


%% Experimental results
load('results/Final_results/2019416-1521_circle_medium_1agent_t45_l216.mat')
N_k = 25;

figure(110)
hold on

[~, sortOrder] = sort(abs(kir));
area_sorted = area(:,sortOrder);

c = getRGBGradient(0.4,0.5,0.6,46);
map = colormap(c);
for i = 1:size(area_sorted,2)
    plot(area_sorted(:,i)/maxArea,'color',c(i,:))
end

bar = colorbar
bar.Ruler.Scale = 'log';
bar.Ruler.MinorTick = 'on';
caxis([0.1 10])

figure(109)
hold on

[meanArea,binKir] = makeMean(kir,N_k,area);

c = getRGBGradient(0.4,0.5,0.6,N_k);
for i = 1:size(meanArea,2)
    plot(meanArea(:,i)/maxArea,'color',c(i,:))
end
bar = colorbar
bar.Ruler.Scale = 'log';
bar.Ruler.MinorTick = 'on';
caxis([0.1 10])

name = 'Experimental results';
title(name)

figure(112)
hold on
c = chir2color(binKir(9:end))
%semilogx(binKir, meanArea(index,:)/maxArea,'k.','markersize',20)

for i = 9:length(binKir(1:end))
    semilogx(binKir(i), meanArea(index,i)/maxArea,'.','color',c(i-8,:),'markersize',20)
end
title('With mean over chirality bins against center of bin')


%%

figure(112)

axis([0.1 10 0 1])
title('Av agenter uppsökt area normerat mot maximal upptäckt area', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 30)

ylabel('Effektivitet', 'Interpreter', 'latex', 'fontsize', 30)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', 30)

legend('Simulering','Exprimentel')



