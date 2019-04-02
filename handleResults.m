%% ploting the results with a loop so we can change parameters

l = 5; % Corresponds to ~5*v*dT for agents in experiments
dT = 1/25;

expName = 'circle_1agent';         %Change name for each new set of data

sourceFile = textscan(fopen(['results/Lab/' expName 'SourceFiles.txt']), '%s','delimiter','\n')
n =size(sourceFile{1},1)

%% HOMOGENOUS
% loop through n XML files

kir = zeros(1,n)
normA = zeros(1,n)
v = zeros(1,n)
for i = 1:n

       file = sourceFile{1}{i}; %we should use same expName here
       [pos_a,~,times] = cut(file,1);
       [kir, v] = getCirality(pos_a,dT,1);
       %spirKir = getChiralitySpiral(pos_a,dT,1,20);
       [squares,normA(i)] = calcArea(pos_a,v,dT,l);
       
       result = [kir(i), normA(i)]              
end


%% COMLEX 2 - with indices 

file = ['results/Lab/' expName 'indices.txt']
indice = load(file);

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
    else j=j+1;
    end 
end 
%% 
%if all works

agent=1;


kir = zeros(1,n)
normA = zeros(1,n)
v = zeros(1,n)
square = zeros(100,n)

for i = 1:n % loop through n XML files
       
       file =  sourceFile{1}{i};
       [pos_a,~,times] = cut(file,1);
       r = zeros(cuts(i),2,size(pos_a,3))

       for j=1:cuts(i)
              r(j,:,1:(indice(new_indice(i,1)+j-1,2)-indice(new_indice(i,1)+j-1,1))+1) = pos_a(agent,:,indice(new_indice(i,1)+j-1,1):indice(new_indice(i,1)+j-1,2)); %picks out cut j from pos_a and makes it agent j in r
              %spirKir(k) = getChiralitySpiral(r,dT,1,20);
       end
       totalTime = size(pos_a,3)*dT;
       
       [kir(i),v(i)] = getComplexCirality(r,dT,1);
       [kir2(i),D_r(i) ,v2(i)] = getKompSpiral(r,dT,1,6,60)
       [square(:,i),normA] = calcArea(pos_a(:,:,250:end),v(i),dT,l,100);
end 
   
%% load result
clear all, close all 

expName = 'Tefat_c1agent';
name= join(['results/Lab/' expName '.txt']);
c= load(name);

kir = c(:,1);
v = c(:,4);
normA = c(:,2);
totalTime = c(:,3);
l = c(:,5);
% D_r = c(:,6);
%% Plot result
kir
[kir, I] = sort(abs(kir));
normA1 = square(I);
k=3;
kir_mm = movmean(kir,k);
normA_mm = movmean(normA1,k);

[kir2, I] = sort(abs(kir2));
normA = normA(I);
kir_mm2 = movmean(kir2,k);
normA_mm2 = movmean(normA,k);
% 
% size_kir=size(kir,1);
% floor_size=floor(size_kir/k)*k;
% 
% normA_m=zeros(floor_size/k,1);
% kir_m=zeros(floor_size/k,1);
% j=1;
% for i=1:size_kir
%     if mod(i,k)==0
%         for m=1:k
%             v(m)=kir(i-m+1);
%             a(m)=normA(i-m+1);
%         end 
%         kir_m(j)=mean(v); 
%         normA_m(j)=mean(a);
%         j=j+1;
%     end 
% end
% clear v
% 
% if size_kir>floor_size
%     rest=size_kir-floor_size
%     for i=1:rest
%         v(i)=kir(floor_size+i-1)
%         a(i)=normA(floor_size+i-1)
%         kir_m(floor_size/k+1)=mean(v);
%         normA_m(floor_size/k+1)=mean(a);
% 
%     end
% end 
%             

figure
semilogx(abs(kir),normA1,'o')
title('no mean')
axis([0.01 10 0 1.1])
% figure
% semilogx(abs(kir_m), normA_m, 'o')
% title('with mean')
% axis([0.01 10 0 1.1])
figure
semilogx(abs(kir), normA_mm, 'o')
title('with moving mean')
axis([0.01 10 0 1.1])


figure
semilogx(abs(kir2), normA_mm2, 'o')
title('with moving mean spiral kirality')
axis([0.01 10 0 1.1])
%% if we want to save the new results
file1 = ['results/Lab/' expName '.txt']; % Name of dataFile

results = zeros(n,6)
results(:,1) = kir;
results(:,2) = normA;
results(:,3) = totalTime;
results(:,4) = v;
results(:,5) = l;
results(:,6) = D_r;

dlmwrite(file1,results);

%%
figure
hold on

[kir2, sortOrder] = sort(abs(kir2))
square = square(:,sortOrder)

color = jet(size(square,2))
for i = 1:size(square,2)
    plot(square(:,i),'color',color(i,:))
end

%%
figure
hold on

sumSquare  = zeros(size(square,1),5)
for k = 1:5
    for j = 1:7
       sumSquare(:,k) = square(:,j+7*(k-1)) + sumSquare(:,k) 
    end
end
sumSquare = sumSquare/7;

color = jet(size(sumSquare,2))
for i = 1:size(sumSquare,2)
    plot(sumSquare(:,i),'color',color(i,:))
end
