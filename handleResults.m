%% ploting the results with a loop so we can change parameters

R = 250;
l = R/10;       % Corresponds to ~5*v*dT for agents in experiments
dT = 1/25;
preTime = 10;   % tid i sekunder innan areaberäkningen börjar
totalTime = 40; % total tid areanberäkningen ska köra efter att den börjat
N = 100;        % antalet tidssteg som calcArea ger tillbaka uppsökt area på
k = 3;          % Hur många movmean medelvärdesbildar på
N_k = 25;     % Antalet bins vi delar upp kiraliteten i 
T = 40;         % Plotta upptäckt area som funktion av kiralitet vid tvärsnitet tiden lika med T s efter pretime

maxArea = pi*R^2
expName = 'circle_medium_1agent';         %Change name for each new set of data

sourceFile = textscan(fopen(['results/Lab/' expName 'SourceFiles.txt']), '%s','delimiter','\n');
n =size(sourceFile{1},1);

%% HOMOGENOUS
% loop through n XML files

kir = zeros(1,n);
normA = zeros(1,n);
v = zeros(1,n);
D_r = zeros(1,n);

for i = 1:n
       file = sourceFile{1}{i}; 
       [pos_a,~,times] = cut(file,1);
       [kir(i),D_r(i) ,v(i)] = getKompSpiral(r,dT,1,6,60);
       [~,normA(i)] = calcArea(pos_a,v(i),dT,l);
       
       result = [kir(i), normA(i)];          
end


%% COMLEX 2 - with indices 
file = ['results/Lab/' expName 'indices.txt']
indice = load(file);

cuts = zeros(1,sum(isnan(indice(:,1))));
new_indice = zeros(sum(isnan(indice(:,1))),2);

film=0;
j=0;
start=1;
for i=1:size(indice,1)
    if isnan(indice(i))
        film=film+1;
        cuts(film)=j; %cuts(i) should be number of cuts in XML file number i 
        new_indice(film,:) = [start start+j-1] ;
        start=start+j+1;
        j=0;
    else
        j=j+1;
    end 
end 
%% 
agent=1;

kir = zeros(1,n);
v = zeros(1,n);
D_r = zeros(1,n);
area = zeros(N,n);

startIndex  =  floor(preTime/dT);
endIndex    =floor((preTime+totalTime)/dT);

for i = 1:n % loop through n XML files
       i
       file =  sourceFile{1}{i};
       [pos_a,~,times] = cut(file,1);
       r = zeros(cuts(i),2,size(pos_a,3)+1); %so as to always have at least one zero 

       for j=1:cuts(i)
              r(j,:,1:(indice(new_indice(i,1)+j-1,2)-indice(new_indice(i,1)+j-1,1))+1) = pos_a(agent,:,indice(new_indice(i,1)+j-1,1):indice(new_indice(i,1)+j-1,2)); %picks out cut j from pos_a and makes it agent j in r
              %spirKir(k) = getChiralitySpiral(r,dT,1,20);
       end
       
       %[kir(i),v(i)] = getComplexCirality(r,dT,1);
       %[kir(i),D_r(i) ,v(i)] = getKompSpiral(r,dT,1,6,60);

       [area(:,i),~] = calcArea(pos_a(:,:,startIndex:endIndex),v(i),dT,l,N);
end 

%% load result
% clear all, close all 

name = join(['results/Lab/' expName '.txt']);
c = load(name);

kir = c(:,1);
% v = c(:,4);
% normA = c(:,2);
% totalTime = c(:,3);
% l = c(:,5);
% D_r = c(:,6);


%% Plotting the area over time for every film
figure
hold on

[~, sortOrder] = sort(abs(kir));
area_sorted = area(:,sortOrder);

color = jet(size(area_sorted,2));
for i = 1:size(area_sorted,2)
    plot(area_sorted(:,i)/maxArea,'color',color(i,:))
end

%% Plotting the mean area of time  over N_k different bins of chirality
figure
hold on

Mi = min(log10(abs(kir)));
Ma = max(log10(abs(kir)));
L = (Ma-Mi)/N_k;

meanArea = zeros(N,N_k);
sumKir = zeros(1,N_k);
binKir = zeros(1,N_k);
count = zeros(1,N_k);
for i = 1:N_k
    for j = 1:size(kir,1)
        if ((i-1)*L + Mi <= log10(abs(kir(j))) && log10(abs(kir(j))) <= Mi + (i)*L )
            meanArea(:,i) = meanArea(:,i) + area(:,j);
            sumKir(i) = sumKir(i) +abs(kir(j));
            count(i) = count(i) +1;
        end
    end
    meanArea(:,i) = meanArea(:,i)/count(i);
    sumKir(i) = sumKir(i)/count(i);
    binKir(i) = 10^(Mi + L*(i+1/2));
end

color = jet(size(meanArea,2));
for i = 1:size(meanArea,2)
    plot(meanArea(:,i)/maxArea,'color',color(i,:))
end


%% Plot the result at end time
index = floor(N*T/totalTime);

[kir_sorted, sortOrder] = sort(abs(kir));
normA1 = area(index,sortOrder);
kir_mm = movmean(kir_sorted,k);
normA_mm = movmean(normA1,k);

figure
semilogx(abs(kir_sorted),normA1./maxArea,'o')
title('no mean')
%axis([0.01 10 0 1.1])

figure
semilogx(abs(kir_sorted), normA_mm./maxArea, 'o')
title('with moving mean on area')
%axis([0.01 10 0 1.1])

figure
semilogx(sumKir, meanArea(index,:)/maxArea,'o')
title('With mean over chirality bins and against mean of chirality')

figure
semilogx(binKir, meanArea(index,:)/maxArea,'o')
title('With mean over chirality bins against center of bin')
%axis([0.01 10 0 0.7])
%% if we want to save the new results
file1 = ['results/Lab/' expName '.txt']; % Name of dataFile

results = zeros(n,6);
results(:,1) = kir;
%results(:,2) = normA;
% results(:,3) = totalTime;
results(:,3) = v;
% results(:,5) = l;
results(:,2) = D_r;

dlmwrite(file1,results);


