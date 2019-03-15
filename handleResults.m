%% ploting the results with a loop so we can change parameters

l = 15; % Corresponds to ~5*v*dT for agents in experiments
dT = 1/25;
n=56;
%complex = false(1,1); % we could merge homogenous and complex

expName = 'hmleft_1agent'; %Change name for each new set of data

%% HOMOGENOUS
% loop through n XML files
for i = 1:n
%      for j = 1:5;           % will not be needed if we name XML files properly
%         k = k+1;
%        if i<10
%            str = join(['0', num2str(j+i*10)]) 
%        else
%           str = num2str(j+i*10)
%        end
       
       str = num2str(i)
       file =  join(['XMLfiles/HomogenLeft_1agent/', str, '_Tracks.xml']) %we should use same expName here
       [pos_a,~,times] = cut(file,1);
       [kir, v] = getCirality(pos_a,dT,1);
       %spirKir = getChiralitySpiral(pos_a,dT,1,20);
       [squares,normA] = calcArea(pos_a,v,dT,l);
       
       result = [kir, normA]
       clc;
        % Sparar:
            % Filnamn/path (i en egen fil)
            % Kiralitet
            % NormArea/tid
            % Total tid
            % hastighet
            % l (size of area elements)
            
        file1 = ['results/Lab/' expName '.txt']; % Name of dataFile
        file2 = ['results/Lab/' expName 'SourceFiles.txt']; % Name of file containing names of XML files
        
        data = [result size(pos_a,3)*dT v l]

        dlmwrite(file1,data,'-append');
        
        fileID = fopen(file2,'a');
        fprintf(fileID,'%-40s\n',file);
        fclose(fileID);
       
end
%% COMPLEX
% loop through n XML files
for i = 1:n
    
       str = num2str(i) %if XMLfiles are named properly
       file =  join(['XMLfiles/HomogenLeft_1agent/', str, '_Tracks.xml'])
       [pos_a,~,times] = cut(file,1);
       
       [r, indice] = splitPositionData(pos_a);
       [kir, v] = getComplexCirality(r,dT,1);
       %spirKir = getChiralitySpiral(pos_a,dT,1,20);
       [squares,normA] = calcArea(pos_a,v,dT,l);
       
       result = [kir, normA]
       
       %now save
       clc;
            
        file1 = ['results/Lab/' expName '.txt']; % Name of dataFile
        file2 = ['results/Lab/' expName 'SourceFiles.txt']; % Name of file containing names of XML files
        
        data = [result size(pos_a,3)*dT v l]

        dlmwrite(file1,data,'-append');

        % complex environment, saves indices where pos_a should be cut for each XML,
        % each XML separated by NaN
        file3 = ['results/Lab/' expName 'indices.txt']
        dlmwrite(file3, indice, '-append')
        
        fileID = fopen(file2,'a');
        fprintf(fileID,'%-40s\n',file);
        fclose(fileID);
       
end

%% load result
clear all, close all 

expName = 'hmleft_1agent';
name= join(['results/Lab/' expName '.txt']);
c= load(name);

[kir, I] = sort(c(:,1));
normA = c((I),2);

k=5;
kir_m = movmean(kir,k);
normA_m = movmean(normA,k);

%% load second result if needed
expName = 'hm1agent';
name= join(['results/Lab/' expName '.txt']);
c= load(name);
[kir2, I] = sort(c(:,1));
normA2 = c((I),2);

k=5; %amount of points for each mean

kir2_m = movmean(kir2,k);
normA2_m = movmean(normA2,k);


%% Plot result

[kir, I] = sort(c(:,1));
normA = c((I),2);

k=5;

kir_m = movmean(kir,k);
normA_m = movmean(normA,k);

figure
semilogx(abs(kir),normA,'o')
title('no mean')
figure
semilogx(abs(kir_m), normA_m, 'o')
title('with mean')

axis([0.01 10 0 1.1])

%% plot joint symmetric result

normA_m_tot=[normA_m; normA2_m];
kir_m_tot = [kir_m; kir2_m];

plot(kir_m_tot, normA_m_tot, 'o')

