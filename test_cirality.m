
file = 'XMLfiles/Homogen_1agent/051_Tracks.xml'
file2 = 'OG.xml'
M = agentTracking(file);
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
clf

pos_a = doublePoint(M,1,6);
plot(squeeze(pos_a(1,1,:)),squeeze(pos_a(1,2,:)),'.')
spirKir = getChiralitySpiral(M,1/25,1,6,30)
%kir = getCirality(M,1/25,1)


%getCirality(pos_a,1/25,1)