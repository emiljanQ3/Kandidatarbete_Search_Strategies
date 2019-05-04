%% Find circle
clf,clc; clear;
% Chose a file with circular trajectory 

load('results/SimulatedTracks/large_w3.30.mat')  % stor
%load('results/SimulatedTracks/medium_w0.187.mat')% mellan
%load('results/SimulatedTracks/small_w0.217.mat') % liten

dT       = 1/25; % Time step
agent    = 1; % which agent/s to plot
digit    = 2;   
agentWidth = 3;
circleWidth = 4;
addRadius = 0.1;

r = R + addRadius; % radius of circle

% Data points on a perfect circle
theta = 0:pi/50:2*pi;
x_val = r * cos(theta);
y_val = r * sin(theta);

% Make two copies of this circle
moveDist = 2*r + r/2;
x_val1 = x_val - moveDist;
x_val3 = x_val + moveDist;

if (R == 1)
    load('results/SimulatedTracks/medium_w0.187.mat')
    M1 = pos_a;
    load('results/SimulatedTracks/medium_w0.9326.mat')
    M2 = pos_a;
    load('results/SimulatedTracks/medium_w4.69.mat')
    M3 = pos_a;
    kir  = abs([0.187; 0.9326; 4.69]); % Medium
elseif(R == 1/1.7)
    load('results/SimulatedTracks/small_w0.217.mat')
    M1 = pos_a;
    load('results/SimulatedTracks/small_w1.55.mat')
    M2 = pos_a;
    load('results/SimulatedTracks/small_w6.27.mat')
    M3 = pos_a;
    kir  = abs([0.217; 1.55; 6.27]); % Liten
else
    load('results/SimulatedTracks/large_w0.207.mat')
    M1 = pos_a;
    load('results/SimulatedTracks/large_w0.5591.mat')
    M2 = pos_a;
    load('results/SimulatedTracks/large_w3.30.mat')
    M3 = pos_a;
    kir = abs([-0.207; 0.5591; -3.30]); % Stor
end
agent = 1;
r = R + addRadius;

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

txt1 = ['$\omega$ = ',num2str(num(1)),' rad/s']; 
txt2 = ['$\omega$ = ',num2str(num(2)),' rad/s'];
txt3 = ['$\omega$ = ',num2str(num(3)),' rad/s'];
     
hold on
p1 = plot(x1,y1,'color',c1,'DisplayName',txt1,'LineWidth',agentWidth)
p2 = plot(x2,y2,'color',c2,'DisplayName',txt2,'LineWidth',agentWidth)
p3 = plot(x3,y3,'color',c3,'DisplayName',txt3,'LineWidth',agentWidth)

plot(x_val1,y_val,'k','LineWidth',circleWidth)
plot(x_val,y_val,'k','LineWidth',circleWidth)
plot(x_val3,y_val,'k','LineWidth',circleWidth)
%plot(x_val,y_val1,'k','LineWidth',6)
%plot(x_val,y_val,'k','LineWidth',4)
%plot(x_val,y_val3,'k','LineWidth',6)

set(gca,'visible','off')
axis equal
xlim([-moveDist-r +moveDist+3*r])
ylim([-1.5*r +1.5*r])

%legend(p1,'Location','north','Interpreter','latex')

%hlegend = legend([p1;p2;p3],'Location','north','Interpreter','latex','FontSize',30)
%hlegend.NumColumns = 3;
%hlegend1 = legend(p1,'Location','west','Interpreter','latex','FontSize',20)
%hlegend2 = legend(p2,'Location','north','Interpreter','latex','FontSize',20)
%hlegend3 = legend(p3,'Location','east','Interpreter','latex','FontSize',20)
%legend(txt1,txt2,txt3)

