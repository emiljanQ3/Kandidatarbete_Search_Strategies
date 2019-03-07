function [kir,normA] = readXML(file,agent)
dT = 1/25;
M = cut(file,agent);

[kir,v] = getCirality(M,dT);
ans = v*dT;

l = 30; % Corresponds to 
[squares,normA] = calcArea(M,v,dT,l);


end