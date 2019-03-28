
file = 'XMLfiles/Homogen_1agent/053_Tracks.xml'
file2 = 'OG.xml'
M = cut(file,1);
%M = cut(file)

figure(1)
%kir = getChiralitySpiral(M,1/25,1,20)


plot(squeeze(M(1,1,:)),squeeze(M(1,2,:)),'.')
axis equal
figure(9111)
pos = doublePoint(M,1,6);
plot(squeeze(pos(1,1,:)),squeeze(pos(1,2,:)),'.')
axis equal
[o v] = getCirality(pos,0.04,1)
%%


pos_a = doublePoint(M,1,6);
%plot(squeeze(pos_a(1,1,:)),squeeze(pos_a(1,2,:)),'.')

[spirKir, D_r ,v] = getChiralitySpiral(M,1/25,1,6,50)

%%


for i = 1:13
    spirKir(1,i) = 10*i;
    [spirKir(2,i), D_r(i) ,v(i)] = getChiralitySpiral(M,1/25,1,6,10*i);
end 
spirKir

w = getCirality(pos_a,0.04,1)
%plot(spirKir,






