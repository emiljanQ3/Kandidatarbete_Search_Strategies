%% Simulated results

load('results/Final_results/2019417-1057_circle_R1_t180_l0156.mat')
    T       = 50;         % vid vilken tidpunkt plottar vi resultatet mot kiralitet
    T2      = 180;         % Hur lång tid plottar vi arean 
    index2  = floor(numAreaDP*T2/measurmentTime);           
    index   = floor(numAreaDP*T/measurmentTime);

    figure(111)
    hold on
    c = chir2color(abs(w));
    map = colormap(c);
    for i = 1:size(meanAreaCovered,2) 
        plot(1:index2, meanAreaCovered(1:index2,i)./maxArea, 'color', c(i,:))
    end
    bar = colorbar
    bar.Ruler.Scale = 'log';
    bar.Ruler.MinorTick = 'on';
    caxis([0.1 10])

    figure(112)
    hold on
    P1 = plot(w,meanAreaCovered(index,:)/maxArea,'b','Linewidth',3)
%     for i = 1:size(meanAreaCovered,2)-1
%         plot([w(i+1) w(i)],[meanAreaCovered(index,i+1), meanAreaCovered(index,i)]/maxArea,'color',c(i,:),'Linewidth',3)
%     end
    set(gca,'xscale','log')
    axis([0.1 10 0 1])
    title('')
%% Experimental results
load('results/Final_results/2019417-1959_circle_medium_1agent_t50_l216.mat')
clear 'length'
N_k = 15;

figure(110)
hold on

[kir_sorted, sortOrder] = sort(abs(kir));
area_sorted = area(:,sortOrder);

c = chir2color(kir_sorted)
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

c = chir2color(binKir)
map = colormap(c);
for i = 1:size(meanArea,2)
    P(i) = plot((1:size(meanArea,1))/2, meanArea(:,i)/maxArea,'color',c(i,:))
end

bar = colorbar
bar.Ruler.Scale = 'log';
bar.Ruler.MinorTick = 'on';
caxis([0.01 10])

name = 'Experimental results';
title(name)

figure(112)
hold on
c = chir2color(binKir)
P2 = semilogx(binKir,meanArea(index,:)/maxArea,'k.','markersize',30)

% for i = 1:length(binKir)
%     semilogx(binKir(i), meanArea(index,i)/maxArea,'.','color',c(i,:),'markersize',30)
% end
title('With mean over chirality bins against center of bin')


%% set font sizes and titel on figur 112
figure(112)

axis([0.1 10 0 1])
title('Experimentell och simulerad data för en agent uppsökt area', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 30)

ylabel('Effektivitet', 'Interpreter', 'latex', 'fontsize', 50)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', 50)

legend([P1 P2], 'Simulering','Experimentell')
%%  set font sizes and titel on figur 109
figure(109)

axis([0 50 0 1])
title('Experimentell data av normerad uppsökt area över tid', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 30)

ylabel('Area ', 'Interpreter', 'latex', 'fontsize', 50)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', 50)
legend(P,binKir)

%% set font sizes and titel on figur 111



