
file = 'XMLfiles/Homogen_1agent/074_Tracks.xml'
file2 = 'OG.xml'

a = 6;  % Antalet pixlar per centimeter

M = cut(file,1);
pos = doublePoint(M,1,7)/6;
M = M/6;

figure(1)
%kir = getChiralitySpiral(M,1/25,1,20)



plot(squeeze(M(1,1,:)),squeeze(M(1,2,:)),'.','Markersize',25)
axis equal

axis([67 120  15 55])
title('Obehandlad positionsdata', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 35)

ylabel('Sträcka (cm)', 'Interpreter', 'latex', 'fontsize', 30)
xlabel('Sträcka (cm)', 'Interpreter', 'latex', 'fontsize', 30)
figure(9111)
plot(squeeze(pos(1,1,:)),squeeze(pos(1,2,:)),'.','Markersize',25)
axis equal
[o v] = getCirality(pos,0.04,1)

axis([67 120  15 55])
title('Interpolerad och filtrerad positionsdata', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 35)

ylabel('Sträcka (cm)', 'Interpreter', 'latex', 'fontsize', 30)
xlabel('Sträcka (cm)', 'Interpreter', 'latex', 'fontsize', 30)

%%


%pos_a = doublePoint(M,1,6);
%plot(squeeze(pos_a(1,1,:)),squeeze(pos_a(1,2,:)),'.')

%[spirKir, D_r ,v] = getChiralitySpiral(M,1/25,1,6,5)

r = splitPositionData(M)

[spirKir2, D_r ,v] = getKompSpiral(r,1/25,1,6,60)
%%


for i = 1:13
    spirKir(1,i) = 10*i;
    [spirKir(2,i), D_r(i) ,v(i)] = getChiralitySpiral(M,1/25,1,6,10*i);
end 
spirKir

w = getCirality(pos_a,0.04,1)
%plot(spirKir,






