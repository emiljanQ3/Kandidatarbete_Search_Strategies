%% Find circle
clc; clear;
% Chose a file with circular trajectory 
file    = 'XMLfiles/Circle_medium_1agent/Medium (1)_Tracks.xml'; 
agent   = 1; 

[pos_a,length,times1] = cut(file,agent); 

maxPos = max(max(pos_a,[],1),[],3);
minPos = min(min(pos_a,[],1),[],3);

sizes = maxPos - minPos;
r = max(sizes)/2 + 7; % radius of circle
x_centre = minPos(1) + max(sizes)/2;
y_centre = minPos(2) + max(sizes)/2;

% Data points on a perfect circle
theta = 0:pi/50:2*pi;
x_val = r * cos(theta) + x_centre;
y_val = r * sin(theta) + y_centre;

% Make two copies of this circle
moveDist = 2*r + r/2;
x_val1 = x_val - moveDist;
x_val3 = x_val + moveDist;

% figure(2)
% hold on
% plot(x_val,y_val,'k','LineWidth',6)
% plot(x_val1,y_val,'k','LineWidth',6)
% plot(x_val3,y_val,'k','LineWidth',6)
% set(gca,'visible','off')
% axis equal

%% Plot trajectrory
clf;

file1    = 'XMLfiles/Circle_medium_1agent/Medium (7)_Tracks.xml'; 
file2    = 'XMLfiles/Circle_medium_1agent/Medium (17)_Tracks.xml'; 
file3    = 'XMLfiles/Circle_medium_1agent/Medium (36)_Tracks.xml'; 
kir      = abs([-0.017954; -0.99264; -4.6946]) 
dT       = 1/25; % Time step
agent    = 1; % which agent/s to plot
digit    = 2;

[M1,length1,times1] = cut(file1,agent); % Turns file into a position matrix, without NaN:s
[M2,length2,times2] = cut(file2,agent); % Turns file into a position matrix, without NaN:s
[M3,length3,times3] = cut(file3,agent); % Turns file into a position matrix, without NaN:s

% Calculate chirality
% myCircle = [pos_a(1,1,200:size(pos_a,3)), pos_a(1,2,200:size(pos_a,3))];
% [r, indice] = splitPositionDataPartitioned(M,700, myCircle);
% [kir,D_r,v] = getKompSpiral(r,dT,1,6,60);

c1 = chir2color(abs(kir(1)));
c2 = chir2color(abs(kir(2)));
c3 = chir2color(abs(kir(3)));

x1 = M1(agent,1,:);
y1 = M1(agent,2,:);
x1 = x1(:,:)'-moveDist;
y1 = y1(:,:)';

x2 = M2(agent,1,:);
y2 = M2(agent,2,:);
x2 = x2(:,:)';
y2 = y2(:,:)';

x3 = M3(agent,1,:);
y3 = M3(agent,2,:);
x3 = x3(:,:)'+moveDist;
y3 = y3(:,:)';

num = round(kir,digit,'significant');


txt1 = ['$\omega$ = ',num2str(num(1)),' rad/s    ']; 
txt2 = ['$\omega$ = ',num2str(num(2)),' rad/s    '];
txt3 = ['$\omega$ = ',num2str(num(3)),' rad/s'];
     

figure(1)
hold on
p1 = plot(x1,y1,'color',c1,'DisplayName',txt1,'LineWidth',3)
p2 = plot(x2,y2,'color',c2,'DisplayName',txt2,'LineWidth',3)
p3 = plot(x3,y3,'color',c3,'DisplayName',txt3,'LineWidth',3)
plot(x_val1,y_val,'k','LineWidth',6)
plot(x_val,y_val,'k','LineWidth',6)
plot(x_val3,y_val,'k','LineWidth',6)
set(gca,'visible','off')
axis equal
xlim([x_centre-moveDist-r x_centre+moveDist+r])
ylim([y_centre-r y_centre+r+100])
hlegend = legend([p1,p2,p3],'Location','north','Interpreter','latex','FontSize',30)
hlegend.NumColumns = 3;