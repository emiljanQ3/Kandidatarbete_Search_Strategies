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
       [kir(k),v(k)] = getCirality2(pos_a,dT);
       [squares,normA(k)] = calcArea(pos_a,v(k),dT,l);

    end
end
figure(117)
semilogx((abs(kir)./v),normA,'o')

figure(118)
semilogx((abs(kir)),normA,'o')


