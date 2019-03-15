%% ploting the results with a loop so we can change parameters

l = 15; % Corresponds to ~5*v*dT for agents in experiments
dT = 1/25;
k= 0;
for i = 1:12
    for j = 1:5 
        k = k+1;
       if i<10
           str = join(['0', num2str(j+i*10)]) 
       else
          str = num2str(j+i*10)
       end
       file =  join(['XMLfiles/Homogen_1agent/', str, '_Tracks.xml'])
       [pos_a,~,times] = cut(file,1);
       [kir(k),v(k)] = getCirality(pos_a,dT,1);
       spirKir(k) = getChiralitySpiral(pos_a,dT,1,20);
       [squares,normA(k)] = calcArea(pos_a,v(k),dT,l);

    end
end

figure(118)
semilogx((abs(kir)),normA,'o')

figure(119)
semilogx((abs(spirKir)),normA,'o')

%% Complex environments
indice = load('results/Lab/c1agentindices.txt')

l = 15; % Corresponds to ~5*v*dT for agents in experiments
dT = 1/25; % Time step
agent = 1; % which agent/s we look at
n=10; %antal 'tiotal' XMLfiler som finns for givet experiment

%cuts=zeros(1,size(indice,1)); 


pos_res = zeros(1,2,) %intiera? hur stor?



film=0; 
j=0;

for i=1:size(indice,1)
    if isnan(indice(i))
        film=film+1;
        cuts(film)=j;
        j=0;
    else j=j+1
    end 
end 
        
k= 0;
for i = 1:n 
    for j = 1:5  %loop over films (k=film), 
        k = k+1;
            if i<10
               str = join(['0', num2str(j+i*10)]) 
           else
              str = num2str(j+i*10)
            end
            
           file =  join(['XMLfiles/c1agent/', str, '_Tracks.xml'])  
           [pos_a,~,times] = cut(file,1);       
           
           

           for j=1:cuts(k)
                  r(j,:,1:(indice(j,2)-indice(j,1)+1)) = pos_a(agent,:,I(1):I(2)); %picks out cut j from pos_a and makes it agent j in r

                  [kir(k),v(k)] = getComplexCirality(r,dT,1);
                  %spirKir(k) = getChiralitySpiral(r,dT,1,20);
                  [squares,normA(k)] = calcArea(pos_a,v(k),dT,l);

           end
    end
end



[pos_a,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

[r, indice] = splitPositionData(pos_a);



[kir,v] = getComplexCirality(r,dT,1);

[squares,normA] = calcArea(pos_a,v,dT,l);

result = [kir, normA]

