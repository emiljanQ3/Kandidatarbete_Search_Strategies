%% Find circle
clf,clc; clear;
% Chose a file with circular trajectory 

%load('results/SimulatedTracks/large_w3.30.mat')  % stor
%load('results/SimulatedTracks/medium_w0.187.mat')% mellan
%load('results/SimulatedTracks/small_w0.217.mat') % liten

load('results/SimulatedTracks/MultiAgent/multiAgent_w0.16.mat')
M1 = pos_a;
load('results/SimulatedTracks/Crosses/crosses_w1.2915.mat')
M2 = pos_a;
load('results/SimulatedTracks/Crosses/crosses_w6.mat')
M3 = pos_a;
kir = [0.16 1.3 6];

dT          = 1/25; % Time step
agent       = 1:2; % which agent/s to plot
digit       = 2;   
agentWidth  = 3;
circleWidth = 4;
addRadius   = 0.5;

r = R + addRadius; % radius of circle

%Data points on a perfect circle
theta = 0:pi/50:2*pi;
x_val = r * cos(theta);
y_val = r * sin(theta);

%Make two copies of this circle
moveDist = 2*r+r/4;
x_val1 = x_val - moveDist;
x_val3 = x_val + moveDist;

% if (R == 1)
%     load('results/SimulatedTracks/medium_w0.187.mat')
%     M1 = pos_a;
%     load('results/SimulatedTracks/medium_w0.9326.mat')
%     M2 = pos_a;
%     load('results/SimulatedTracks/medium_w4.69.mat')
%     M3 = pos_a;
%     kir  = abs([0.187; 0.9326; 4.69]); % Medium
% elseif(R == 1/1.7)
%     load('results/SimulatedTracks/small_w0.217.mat')
%     M1 = pos_a;
%     load('results/SimulatedTracks/small_w1.55.mat')
%     M2 = pos_a;
%     load('results/SimulatedTracks/small_w6.27.mat')
%     M3 = pos_a;
%     kir  = abs([0.217; 1.55; 6.27]); % Liten
% else
%     load('results/SimulatedTracks/large_w0.207.mat')
%     M1 = pos_a;
%     load('results/SimulatedTracks/large_w0.5591.mat')
%     M2 = pos_a;
%     load('results/SimulatedTracks/large_w3.30.mat')
%     M3 = pos_a;
%     kir = abs([-0.207; 0.5591; -3.30]); % Stor
% end
agent = 1:2;
r = R + addRadius;

c1 = chir2color(abs(kir(1)));
c2 = chir2color(abs(kir(2)));
c3 = chir2color(abs(kir(3)));

x1 = M1(agent,1,:);
y1 = M1(agent,2,:);
x1 = x1(:,:)'-moveDist;
y1 = y1(:,:)';
agent=1;
x2 = M2(agent,1,:);
y2 = M2(agent,2,:);
x2 = x2(:,:)'+30;
y2 = y2(:,:)'+30;

x3 = M3(agent,1,:);
y3 = M3(agent,2,:);
x3 = x3(:,:)'+moveDist;
y3 = y3(:,:)'+30;

num = round(kir,digit,'significant');

txt1 = ['$\omega$ = ',num2str(num(1)),' rad/s']; 
txt2 = ['$\omega$ = ',num2str(num(2)),' rad/s'];
txt3 = ['$\omega$ = ',num2str(num(3)),' rad/s'];
     
hold on
% maxPos = max(max(M3,[],1),[],3);
% minPos = min(min(M3,[],1),[],3);
% Plot crosses
% maxPos = [210 120];
% minPos = [-150 0];
% maxmax = max(abs([maxPos minPos]));
% plotSize = ceil((maxmax) / L);
% for i = -plotSize():plotSize-1
%     for j = -plotSize:plotSize-1
%         for k = 1:size(obstacle, 3)
%             plot(L*obstacle(:,1,k)+j*L,L*obstacle(:,2,k)+i*L, 'k', 'LineWidth', 1)
%             
%         end
%     end
% end


hold on
p1 = plot(x1,y1,'color',c1,'DisplayName',txt1,'LineWidth',agentWidth)
p2 = plot(x2,y2,'color',c2,'DisplayName',txt2,'LineWidth',agentWidth)
p3 = plot(x3,y3,'color',c3,'DisplayName',txt3,'LineWidth',agentWidth)

plot(x_val1,y_val,'k','LineWidth',circleWidth)
plot(x_val,y_val,'k','LineWidth',circleWidth)
plot(x_val3,y_val,'k','LineWidth',circleWidth)
% plot(x_val,y_val,'k','LineWidth',6)
% plot(x_val,y_val,'k','LineWidth',4)
% plot(x_val,y_val3,'k','LineWidth',6)

set(gca,'visible','off')
axis equal
xlim([min(x_val1) max(x_val3)])
ylim([min(y_val) max(y_val)])

%legend(p1,'Location','north','Interpreter','latex')
% hlegend = legend([p1,p2,p3],'Location','north','Interpreter','latex','FontSize',20)
% hlegend.NumColumns = 3;
a = axes('position',get(gca,'position'),'visible','off')
b = copyobj(a,gcf);
hlegend1 = legend(a,p1,'Location','northwestoutside','Interpreter','latex','FontSize',20)
hlegend2 = legend(b,p2,'Location','north','Interpreter','latex','FontSize',20)
hlegend3 = legend(p3,'Location','northeast','Interpreter','latex','FontSize',20)
%legend(txt1,txt2,txt3)

%% MultiAgents

clf,clc; clear;
% Chose a file with circular trajectory 

%load('results/SimulatedTracks/large_w3.30.mat')  % stor
%load('results/SimulatedTracks/medium_w0.187.mat')% mellan
%load('results/SimulatedTracks/small_w0.217.mat') % liten

load('results/SimulatedTracks/MultiAgent/multiAgent_w1.8.mat');
M1 = pos_a;
M1 = removezeros(M1);
kir1 = w;
load('results/SimulatedTracks/MultiAgent/multiAgent_w-4.6.mat');
M2 = pos_a;
M2 = removezeros(M2);
kir2 = w;
load('results/SimulatedTracks/MultiAgent/multiAgent_w-0.44.mat');
M3 = pos_a;
M3 = removezeros(M3);
kir3 = w;
load('results/SimulatedTracks/MultiAgent/multiAgent_w0.5.mat');
M4 = pos_a;
M4 = removezeros(M4);
kir4 = w;
load('results/SimulatedTracks/MultiAgent/multiAgent_w0.8.mat');
M5 = pos_a;
M5 = removezeros(M5);
kir5 = w;
load('results/SimulatedTracks/MultiAgent/multiAgent_w2.6.mat');
M6 = pos_a;
M6 = removezeros(M6);
kir6 = w;

kir = [kir1' kir2' kir3' kir4' kir5' kir6'];

dT          = 1/25; % Time step
agent       = 1:2; % which agent/s to plot
digit       = 2;   
agentWidth  = 3;
circleWidth = 4;
addRadius   = 0.2;
fontSize    = 12;

r = R + addRadius; % radius of circle

%Data points on a perfect circle
theta = 0:pi/50:2*pi;
x_val = r * cos(theta);
y_val = r * sin(theta);

%Make two copies of this circle
moveDist = 2*r+r/4;
x_val1 = x_val - 3*moveDist;
x_val2 = x_val - 2*moveDist;
x_val3 = x_val - moveDist;
x_val4 = x_val;
x_val5 = x_val + moveDist;
x_val6 = x_val + 2*moveDist;

r = R + addRadius;

c1 = chir2color(abs(kir(1,1)));
c2 = chir2color(abs(kir(2,1)));
c3 = chir2color(abs(kir(1,2)));
c4 = chir2color(abs(kir(2,2)));
c5 = chir2color(abs(kir(1,3)));
c6 = chir2color(abs(kir(2,3)));
c7 = chir2color(abs(kir(1,4)));
c8 = chir2color(abs(kir(2,4)));
c9 = chir2color(abs(kir(1,5)));
c10 = chir2color(abs(kir(2,5)));
c11 = chir2color(abs(kir(1,6)));
c12 = chir2color(abs(kir(2,6)));

[x1,y1] = getXYVectors(M1,-3*moveDist);
[x2,y2] = getXYVectors(M2,-2*moveDist);
[x3,y3] = getXYVectors(M3,-moveDist);
[x4,y4] = getXYVectors(M4,0);
[x5,y5] = getXYVectors(M5,moveDist);
[x6,y6] = getXYVectors(M6,2*moveDist)

y4(:,2) = y4(:,2) + 0.05;

num = round(kir,digit,'significant');

txt1 = ['$\omega$ = ',num2str(num(1,1)),' rad/s$\quad$']; 
txt2 = ['$\omega$ = ',num2str(num(2,1)),' rad/s$\quad$'];
txt3 = ['$\omega$ = ',num2str(num(1,2)),' rad/s$\quad$'];
txt4 = ['$\omega$ = ',num2str(num(2,2)),' rad/s$\quad$'];
txt5 = ['$\omega$ = ',num2str(num(1,3)),' rad/s$\quad$'];
txt6 = ['$\omega$ = ',num2str(num(2,3)),' rad/s$\quad$'];
txt7 = ['$\omega$ = ',num2str(num(1,4)),' rad/s$\quad$']; 
txt8 = ['$\omega$ = ',num2str(num(2,4)),' rad/s$\quad$'];
txt9 = ['$\omega$ = ',num2str(num(1,5)),' rad/s$\quad$'];
txt10 = ['$\omega$ = ',num2str(num(2,5)),' rad/s$\quad$'];
txt11 = ['$\omega$ = ',num2str(num(1,6)),' rad/s$\quad$'];
txt12 = ['$\omega$ = ',num2str(num(2,6)),' rad/s$\quad$'];

hold on
p1 = plot(x1(:,1),y1(:,1),'color',c1,'DisplayName',txt1,'LineWidth',agentWidth);
p2 = plot(x1(:,2),y1(:,2),'color',c2,'DisplayName',txt2,'LineWidth',agentWidth);
p3 = plot(x2(:,1),y2(:,1),'color',c3,'DisplayName',txt3,'LineWidth',agentWidth);
p4 = plot(x2(:,2),y2(:,2),'color',c4,'DisplayName',txt4,'LineWidth',agentWidth);
p5 = plot(x3(:,1),y3(:,1),'color',c5,'DisplayName',txt5,'LineWidth',agentWidth);
p6 = plot(x3(:,2),y3(:,2),'color',c6,'DisplayName',txt6,'LineWidth',agentWidth);
p7 = plot(x4(:,1),y4(:,1),'color',c7,'DisplayName',txt7,'LineWidth',agentWidth);
p8 = plot(x4(:,2),y4(:,2),'color',c8,'DisplayName',txt8,'LineWidth',agentWidth);
p9 = plot(x5(:,1),y5(:,1),'color',c9,'DisplayName',txt9,'LineWidth',agentWidth);
p10 = plot(x5(:,2),y5(:,2),'color',c10,'DisplayName',txt10,'LineWidth',agentWidth);
p11 = plot(x6(:,1),y6(:,1),'color',c11,'DisplayName',txt11,'LineWidth',agentWidth);
p12 = plot(x6(:,2),y6(:,2),'color',c12,'DisplayName',txt12,'LineWidth',agentWidth);

% Plot collisions
[time1] = checkCollision(x1,y1);
[time2] = checkCollision(x2,y2);
[time3] = checkCollision(x3,y3);
[time4] = checkCollision(x4,y4);
[time5] = checkCollision(x5,y5);
[time6] = checkCollision(x6,y6);

plot(x_val1,y_val,'k','LineWidth',circleWidth)
plot(x_val2,y_val,'k','LineWidth',circleWidth)
plot(x_val3,y_val,'k','LineWidth',circleWidth)
plot(x_val4,y_val,'k','LineWidth',circleWidth)
plot(x_val5,y_val,'k','LineWidth',circleWidth)
plot(x_val6,y_val,'k','LineWidth',circleWidth)
% plot(x_val,y_val,'k','LineWidth',6)
% plot(x_val,y_val,'k','LineWidth',4)
% plot(x_val,y_val3,'k','LineWidth',6)

set(gca,'visible','off')
axis equal
xlim([min(x_val1) max(x_val6)])
ylim([min(y_val) max(y_val)])

%legend(p1,'Location','north','Interpreter','latex')
% hlegend = legend([p1,p2,p3],'Location','north','Interpreter','latex','FontSize',20)
% hlegend.NumColumns = 3;
%a = axes('position',get(gca,'position'),'visible','off')
% b = copyobj(a,gcf);
% hlegend1 = legend(a,p1,'Location','northwestoutside','Interpreter','latex','FontSize',20)
% hlegend2 = legend(b,p2,'Location','north','Interpreter','latex','FontSize',20)
% hlegend3 = legend(p3,'Location','northeast','Interpreter','latex','FontSize',20)
%legend(txt1,txt2,txt3)
a = axes('position',get(gca,'position'),'visible','off');
b = copyobj(a,gcf);
c = copyobj(b,gcf);
d = copyobj(c,gcf);
e = copyobj(b,gcf);
f = copyobj(c,gcf);
% hlegend1 = legend(a,[p1,p2],'Location','northwest','Interpreter','latex','FontSize',20);
% hlegend2 = legend(b,[p3,p4],'Location','north','Interpreter','latex','FontSize',20);
% hlegend3 = legend(c,[p5,p6],'Location','northeast','Interpreter','latex','FontSize',20);
hl = legend([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12],'Location','northeast','Interpreter','latex','FontSize',fontSize)
hl.NumColumns = 6;

dist = 0.132;
annotation('textbox',[0.146 0.28 0.1 0.1],'String',time1,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.146+dist 0.28 0.1 0.1],'String',time2,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.146+2*dist-0.001 0.28 0.1 0.1],'String',time3,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.146+3*dist-0.007 0.28 0.1 0.1],'String',time4,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.146+4*dist-0.005 0.28 0.1 0.1],'String',time5,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.146+5*dist-0.005 0.28 0.1 0.1],'String',time6,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);


%%
function [data_out] = removezeros(data_in)
    nonzeroindices = find(nonzeros(data_in(1,1,:)));
    data_out = data_in(:,:,nonzeroindices);
end

%%
function [x,y] = getXYVectors(M,moveDist)
agent = 1:2
x = M(agent,1,:);
y = M(agent,2,:);
x = x(:,:)'+moveDist;
y = y(:,:)';
end

%%
function [time] = checkCollision(x,y)
dT = 0.04;
if (length(x)<1500)
    midP = [(x(end,1)+x(end,2))/2, (y(end,1)+y(end,2))/2];
    plot(midP(1),midP(2),'kx','LineWidth',3,'MarkerSize',20);
    time = [num2str(length(x)*dT) ' sekunder'];
else
    time = 'Maximal tid';
end
end