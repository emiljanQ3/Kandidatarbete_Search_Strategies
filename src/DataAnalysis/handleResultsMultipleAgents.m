%% ploting the results with a loop so we can change parameters
clc
R = 270/2;
l = R/20;       % Corresponds to ~5*v*dT for agents in experiments
dT = 1/25;
preTime = 10;   % tid i sekunder innan areaberäkningen börjar
totalTime = 45; % total tid areanberäkningen ska köra efter att den börjat
N = 100;        % antalet tidssteg som calcArea ger tillbaka uppsökt area på
k = 3;          % Hur många movmean medelvärdesbildar på
N_k = 25;     % Antalet bins vi delar upp kiraliteten i
T = 40;         % Plotta upptäckt area som funktion av kiralitet vid tvärsnitet tiden lika med T s efter pretime
agent = 1:2;

maxArea = pi*R^2;
expName = 'circle_medium_2agent';         % Change name for each new set of data

sourceFile = textscan(fopen(['results/Lab/' expName 'SourceFiles.txt']), '%s','delimiter','\n');
n =size(sourceFile{1},1);

%% COMPLEX 2 - with indices
file = ['results/Lab/' expName 'indices.txt']
indice = load(file);

% Number of sets of indices
cuts = zeros(2,sum(isnan(indice(:,1)))/2);
% Contains indices of start and end position for a film in indice
% First agent at new_indice(film,1:2), second at new_indice(film,3:4)
new_indice = zeros(sum(isnan(indice(:,1)))/2,4);

film=0;
j=0;
start=1;
bool = false; % False if not both agents indices are found
% Finds all values for new_indice
for i=1:size(indice,1)
    if (isnan(indice(i)) && bool == false) % End of one set
        film=film+1;
        cuts(1,film)=j; %cuts(i) should be number of cuts in XML file number i
        new_indice(film,1:2) = [start start+j-1] ;
        start=start+j+1;
        j=0;
        bool = true;
    elseif (isnan(indice(i)) && bool == true)
        cuts(2,film)=j;
        new_indice(film,3:4) = [start start+j-1] ;
        start=start+j+1;
        j=0;
        bool = false;
    else
        % j is number of cuts in one video
        j=j+1;
    end
end


kir1 = zeros(1,n);
v1 = zeros(1,n);
D_r1 = zeros(1,n);
kir2 = zeros(1,n);
v2 = zeros(1,n);
D_r2 = zeros(1,n);
% area = zeros(N,n);

startIndex  =  floor(preTime/dT); % First frame
endIndex    = floor((preTime+totalTime)/dT); % Last frame

% Saves the total time of each video
total_time = zeros(n,1)

for i = 1:n % loop through n XML files
    file =  sourceFile{1}{i};
    [pos_a,~,times] = cut(file,agent);
    total_time(i,1) = times(1,2); % Assumes both agents exist for the same amount of time
    pos1 = zeros(cuts(1,i),2,size(pos_a,3)+1); % so as to always have at least one zero
    pos2 = zeros(cuts(2,i),2,size(pos_a,3)+1);
    for j=1:cuts(1,i)
        %picks out cut j from pos_a and makes it agent j in r
        pos1(j,:,1:(indice(new_indice(i,1)+j-1,2)-indice(new_indice(i,1)+j-1,1))+1) = pos_a(1,:,indice(new_indice(i,1)+j-1,1):indice(new_indice(i,1)+j-1,2));
    end
    
    for j=1:cuts(2,i)
        %picks out cut j from pos_a and makes it agent j in r
        pos2(j,:,1:(indice(new_indice(i,3)+j-1,2)-indice(new_indice(i,3)+j-1,1))+1) = pos_a(2,:,indice(new_indice(i,3)+j-1,1):indice(new_indice(i,3)+j-1,2));
    end
    
    
    
    %[kir(i),v(i)] = getComplexCirality(r,dT,1);
    
    %w = waitforbuttonpress;
    
    % Calculate chirality, D_r and speed separately for each agent
    [kir1(i),D_r1(i) ,v1(i)] = getKompSpiral(pos1,dT,1,6,60);
    [kir2(i),D_r2(i) ,v2(i)] = getKompSpiral(pos2,dT,1,6,60);
    %[kir1(i) D_r(i) v(i)]
    
end

%% load result
% clear all, close all
% Not to be used later, I think

name = join(['results/Lab/' expName '.txt']);
c = load(name);

kir_saved = c(:,1);
D_r_saved = c(:,2);
v_saved = c(:,3);
time_saved = c(:,4);
% totalTime = c(:,3);
% l = c(:,5);

for i = 1:size(kir_saved)/2
   kir1(i) = kir_saved(2*i-1);
   kir2(i) = kir_saved(2*i);
   time(i) = time_saved(2*i-1);
end
%(kir1,kir2,time,'.')

%Mirror points
time = [time, time];
kir1 = [kir1, kir2];
kir2 = [kir2, kir1];

%% Plot datapoints binned
% Scatters data points, time/efficiency is colour coded
kirRange = 5;
ax_Font = 40;
gca_Font = 30;

[time_bin, kir1_bin, kir2_bin] = linearBin2(kir1, kir2, 40, time, [-5 5]);

k = 1;
for i = 1:length(kir1_bin)
    for j = 1:length(kir2_bin)
        if ~isnan(time_bin(i,j))
            time_bin2(k) = time_bin(i,j);
            kir1_bin2(k) = kir1_bin(i);
            kir2_bin2(k) = kir2_bin(j);
            
            k = k + 1;
        end
    end
end

hold on
[total_time_sorted, sortOrder] = sort(time_bin2);
kir1_sorted = kir1_bin2(sortOrder);
kir2_sorted = kir2_bin2(sortOrder);

% %plot 1
% figure
% hold on
% img = imread('Images/svartvitTid.png');
% image('CData', img, 'XData', [kirRange ,-kirRange], 'YData', [-kirRange, kirRange])
% scatter(kir1_sorted,kir2_sorted, [],total_time_sorted)
% 
% %Formatting 1
% axis('square')
% axis([-kirRange kirRange -kirRange kirRange])


%Plot 2
efficiency = 1./total_time_sorted;
figure
hold on
img = imread('Images/svartvitEffektivitet.png');
image('CData', img, 'XData', [kirRange ,-kirRange], 'YData', [-kirRange, kirRange])
scatter(kir1_sorted,kir2_sorted, 'k.', 'SizeData', 1200)
scatter(kir1_sorted,kir2_sorted, [],efficiency, '.', 'SizeData', 1000)


%Formatting 2
axis('square')
axis([-kirRange kirRange -kirRange kirRange])
set(gca, 'fontsize', gca_Font)
xticks(-kirRange:1:kirRange)
yticks(-kirRange:1:kirRange)
set(gca,'Layer','top','GridAlpha', 0.0)

ylabel('Kiralitet agent A (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet agent B (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
bar = colorbar;
set(get(bar,'label'),'string','Effektivitet (s$^{-1}$)', 'Interpreter', 'latex', 'fontsize', ax_Font);


%% Plot datapoints
% Scatters data points, time/efficiency is colour coded
kirRange = 3;

hold on
[total_time_sorted, sortOrder] = sort(time);
kir1_sorted = kir1(sortOrder);
kir2_sorted = kir2(sortOrder);

%plot 1
figure(21)
hold on
img = imread('Images/svartvitTid.png');
image('CData', img, 'XData', [kirRange ,-kirRange], 'YData', [-kirRange, kirRange])
scatter(kir1_sorted,kir2_sorted, [],total_time_sorted)

%Formatting 1
axis('square')
axis([-kirRange kirRange -kirRange kirRange])

%Plot 2
efficiency = 1./total_time_sorted;
figure(22)
hold on
img = imread('Images/svartvitEffektivitet.png');
image('CData', img, 'XData', [kirRange ,-kirRange], 'YData', [-kirRange, kirRange])
scatter(kir1_sorted,kir2_sorted, 'k.', 'SizeData', 1200)
scatter(kir1_sorted,kir2_sorted, [],efficiency, '.', 'SizeData', 1000)
%plot([-3, 3], [-1,-1], 'r', 'Linewidth', 2)

%Formatting 2
axis('square')
axis([-kirRange kirRange -kirRange kirRange])
set(gca, 'fontsize', gca_Font)
xticks(-kirRange:1:kirRange)
yticks(-kirRange:1:kirRange)
set(gca,'Layer','top','GridAlpha', 0.0)

ylabel('Kiralitet agent A (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet agent B (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
bar = colorbar;
set(get(bar,'label'),'string','Effektivitet (s$^{-1}$)', 'Interpreter', 'latex', 'fontsize', ax_Font);

%% Tvärsnitt
rangeKir1 = [-3 3];
rangeKir2 = [-1.2 -0.7];

j = 1;
for i = 1:length(kir1)
    
    if ~(kir1(i) < rangeKir1(1) || kir1(i) > rangeKir1(2) || kir2(i) < rangeKir2(1) || kir2(i) > rangeKir2(2))
        kirLine1(j) = kir1(i);
        kirLine2(j) = kir2(i);
        timeLine(j) = time(i);
        
        j = j + 1;
        continue;
    end
    %time(i) = NaN;
    %kir2(i) = NaN;
    %kir1(i) = NaN;
end
 
[binnedTimes, binKir1] = linearBin(kirLine1, 18, timeLine, rangeKir1);
effic = 1./binnedTimes;

figure(112)
hold on
c = chir2color(abs(binKir1));
map = colormap(chir2color(10.^(-1:0.02:1)));
P2 = plot(binKir1,effic,'k.','markersize',50)

for i = 1:length(binKir1)
    plot(binKir1(i), effic(i),'.','color',c(i,:),'markersize',40)
end
bar = colorbar('Location', 'northoutside', 'Ticks', ([0.1 4.95 10]), 'TickLabels',{'0.1', '1','10'});
bar.Ruler.MinorTick = 'off';
caxis([0.1 10])

%%
figure(1338)
hold on
effic2 = 1./timeLine;
scatter(kirLine1, effic2);


%% Plot interpolated data as surface plot - not done...
xlin = linspace(min(kir1),max(kir1),50);
ylin = linspace(min(kir2),max(kir2),50);
[X,Y] = meshgrid(xlin,ylin);

f = scatteredInterpolant(kir1,kir2,total_time);
Z = f(X,Y);

%% if we want to save the new results
% Not done for multiple agents
file1 = ['results/Lab/' expName '.txt']; % Name of dataFile

results = zeros(n,6);
results(:,1) = kir1;
%results(:,2) = normA;
% results(:,3) = totalTime;
results(:,3) = v;
% results(:,5) = l;
results(:,2) = D_r;

dlmwrite(file1,results);

%% Run to save workspace
dateTime = clock;
R_s = num2str(R);
l_s = num2str(l);
time_s = num2str(totalTime);

filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''),'_',expName, '_t', time_s, '_l', l_s([1:2,4:end]));
path = strcat(pwd, '/results/Final_results/', filename)
save(path)


