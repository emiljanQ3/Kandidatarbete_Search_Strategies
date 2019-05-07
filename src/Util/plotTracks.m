%% Find circle
clc; clear;
% Chose a file with circular trajectory 
path = 'XMLfiles/Circle_medium_1agent/Medium '
file    = [path '(4)_Tracks.xml']; 
agent   = 1; 

[pos_a,length,times1] = cut(file,agent); 

maxPos = max(max(pos_a,[],1),[],3);
minPos = min(min(pos_a,[],1),[],3);

sizes = maxPos - minPos;
r = max(sizes)/2+20; % radius of circle
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

% y_val1 = y_val - moveDist;
% y_val3 = y_val + moveDist;

% figure(2)
% hold on
% plot(x_val,y_val,'k','LineWidth',6)
% plot(x_val1,y_val,'k','LineWidth',6)
% plot(x_val3,y_val,'k','LineWidth',6)
% set(gca,'visible','off')
% axis equal

%% Plot trajectrory
clf;
% Liten: 101, 77, 81
% Medium: 4, 17, 36
% Stor: 17, 12, 35

file1    = file; 
file2    = [path '(17)_Tracks.xml']; 
file3    = [path '(36)_Tracks.xml']; 
%kir      = abs([-0.20704; 0.49383; -3.3031]) % Stor
kir      = abs([0.18783; -0.99264; -4.6946]) % Medium
%kir      = abs([-0.2172; -1.3443; -6.272]) % Liten
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

x1 = M1(agent,1,100:end);
y1 = M1(agent,2,100:end);
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
     
img = imread('hmMedBack.PNG');
img_size = r+20;
lightGrey1   = [0.85 0.85 0.85];

figure(1)
hold on
image('CData',img,'XData',[x_centre-img_size-moveDist x_centre+img_size-moveDist],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size x_centre+img_size],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size+moveDist x_centre+img_size+moveDist],'YData',[y_centre-img_size y_centre+img_size])

p1 = plot(x1,y1,'color',c1,'DisplayName',txt1,'LineWidth',3)
p2 = plot(x2,y2,'color',c2,'DisplayName',txt2,'LineWidth',3)
p3 = plot(x3,y3,'color',c3,'DisplayName',txt3,'LineWidth',3)

plot(x_val1,y_val,'Color',lightGrey1,'LineWidth',4)
plot(x_val,y_val,'Color',lightGrey1,'LineWidth',4)
plot(x_val3,y_val,'Color',lightGrey1,'LineWidth',4)
%plot(x_val,y_val1,'k','LineWidth',6)
%plot(x_val,y_val,'k','LineWidth',4)
%plot(x_val,y_val3,'k','LineWidth',6)

set(gca,'visible','off')
axis equal
xlim([y_centre-moveDist-r y_centre+moveDist+3*r])
ylim([x_centre-1.5*r x_centre+1.5*r])
%hlegend = legend([p1,p2,p3],'Location','north','Interpreter','latex','FontSize',30)
%hlegend.NumColumns = 3;
hlegend1 = legend([p1],'Location','west','Interpreter','latex','FontSize',20)
hlegend2 = legend([p2],'Location','north','Interpreter','latex','FontSize',20)
hlegend3 = legend([p3],'Location','east','Interpreter','latex','FontSize',20)
%legend(txt1,txt2,txt3)