%% Parameters
clc,,clf
digit = 2;              % Antalet signifikanta siffror vi skriver ut f�r kiralitet i legend
N_K = 15;               % Antalet kiralitetsbins vi delar upp datan i           

T       = 90;           % vid vilken tidpunkt plottar vi resultatet mot kiralitet
T2      = 90;           % Hur lång tid plottar vi arean 

index2  = floor(numAreaDP*T2/measurmentTime);
index1   = floor(numAreaDP*T/measurmentTime);

linewidth = 4;          % Hur breda plottade linjer ska vara
markSize = 40;          % Hur stor plottade datapunkter ska vara

ax_Font = 50;           % Fontsizes on axis
tit_Font = 30;          % Fontsizes on titels
%% Experimental results
clc
load('results/Final_results/2019417-118_circle_large_1agent_t50_l216.mat')
clear 'length'


c1 = chir2color(10.^(-1:0.02:1));

% In figure(110) we are plotting al the kiralitys should not be in report
figure(110)
hold on
[kir_sorted, sortOrder] = sort(abs(kir));
area_sorted = area(:,sortOrder);

c = chir2color(kir_sorted);
map = colormap(c1);
for i = 1:size(area_sorted,2)
    plot(area_sorted(:,i)/maxArea,'color',c(i,:))
end
bar = colorbar;
bar.Ruler.Scale = 'log';
bar.Ruler.MinorTick = 'on';
caxis([0.1 10])

% In figure(109) we are plotting the kiralitys which are closest to the
% chirality bins for the given number of chiralitybins goes in report
figure(109)
hold on

[meanArea,binKir] = makeMean(kir,N_K,area);

c = chir2color(binKir);
map = colormap(c1);
format short
for i = 1:size(meanArea,2)
    num = round(binKir(i),digit,'significant');
    txt = ['$\omega$ = ',num2str(num),' rad/s'];
    P(i) = plot((1:size(meanArea,1))*9/10, meanArea(:,i)/maxArea,'color',c(i,:),'DisplayName',txt,'Linewidth', linewidth)
end

bar = colorbar;
bar.Ruler.Scale = 'log';
bar.Ruler.MinorTick = 'on';
caxis([0.1 10]);

% In figure(112) we are plotting area over chiraleties both for sim and exp
% at the time T given in parameters
figure(112)
hold on
c = chir2color(binKir);
map = colormap(c1);
P2 = semilogx(binKir,meanArea(index1,:)/maxArea,'k.','markersize',50)
set(gca,'xscale','log')
for i = 1:length(binKir)
    semilogx(binKir(i), meanArea(index1,i)/maxArea,'.','color',c(i,:),'markersize',markSize)
end
bar = colorbar;
bar.Ruler.Scale = 'log';
bar.Ruler.MinorTick = 'on';
caxis([0.1 10])
%% Simulated results
load('results/Final_results/2019417-1359_circle_R17_t360_l0156.mat')
clear 'min'
clear 'max'
    % In figure(111) we plot al simulated ciralitys not in report
    figure(111)
    hold on
    c = chir2color(abs(w));
    map = colormap(c);
    for i = 1:size(meanAreaCovered,2) 
        plot(1:index2, meanAreaCovered(1:index2,i)./maxArea, 'color', c(i,:))
    end
    bar = colorbar;
    bar.Ruler.Scale = 'log';
    bar.Ruler.MinorTick = 'on';
    caxis([0.1 10])
    
    ind = zeros(length(binKir),1);
    for i=1:length(binKir)
        W = abs(w-binKir(i));
        [w2, ind(i)] = min(W);
    end
    
    % In figure(113) we plot area over time for all the chirality binsgoes
    % in report
    figure(113)
    map = colormap(c);
    hold on
    for i = ind' 
        num = round(w(i),digit,'significant');
        txt = ['$\omega$ = ',num2str(num),' rad/s'];
        plot(1:index2, meanAreaCovered(1:index2,i)./maxArea, 'color', c(i,:),'DisplayName',txt,'Linewidth',linewidth)
    end
    bar = colorbar;
    bar.Ruler.Scale = 'log';
    bar.Ruler.MinorTick = 'on';
    caxis([0.1 10])
    
    figure(112)
    hold on
    P1 = plot(w,meanAreaCovered(index1,:)/maxArea,'k','Linewidth',linewidth)
%     for i = 1:size(meanAreaCovered,2)-1
%         plot([w(i+1) w(i)],[meanAreaCovered(index,i+1), meanAreaCovered(index,i)]/maxArea,'color',c(i,:),'Linewidth',3)
%     end
    set(gca,'xscale','log')
%% set font sizes and titel on figur 112
figure(112)

axis([0.1 10 0 1])
title('Experimentell och simulerad uppsökt area för en agent', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Normerad area', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)

%legend([P1 P2], 'Simulering','Experimentell','Location','northwest')
%%  set font sizes and titel on figur 109
figure(109)

axis([0 T2 0 1])
title('Experimentell data av normerad uppsökt area över tid', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Normerad area ', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Tid (s)', 'Interpreter', 'latex', 'fontsize', ax_Font)

%hlegend = legend('show','Location','northwest','Interpreter','latex')
%hlegend.NumColumns=2;      % g�r legenden i tv� kolumner

%% set font sizes and titel on figur 113
figure(113)

axis([0 T2 0 1])
title('Simulerad data av normerad uppsökt area över tid', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Normerad area', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Tid (s)', 'Interpreter', 'latex', 'fontsize', ax_Font)

%hlegend = legend('show','Location','northwest','Interpreter','latex');
%hlegend.NumColumns = 2;      % gör legenden i 2 kolumner


