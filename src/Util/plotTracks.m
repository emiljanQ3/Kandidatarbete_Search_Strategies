%% Find circle

% For multiagents make more circles
clc; clear;
% exp = 'large';
exp = 'medium';
% exp = 'small';

if strcmp(exp,'large')
    file = 'XMLfiles/Circle_large_1agent/Stor (17)_Tracks.xml';
    r_add = 35;
    x_add = -10;
    y_add = 0;
    move_add = 60;
elseif strcmp(exp,'medium')
    file = 'XMLfiles/Circle_medium_1agent/Medium (4)_Tracks.xml';
    r_add = 17;
    x_add = 13;
    y_add = 10;
    move_add = 55;
elseif strcmp(exp,'small')
    file = 'XMLfiles/Circle_small_1agent/Liten (101)_Tracks.xml';
    r_add = 10;
    x_add = 0;
    y_add = 0;
    move_add = 27;
end
% Chose a file with circular trajectory 
% path = 'XMLfiles/Circle_medium_1agent/Medium '
% file    = [path '(17)_Tracks.xml']; 

agent   = 1; 

[pos_a] = cut(file,agent); 

maxPos = max(max(pos_a,[],1),[],3);
minPos = min(min(pos_a,[],1),[],3);

sizes = maxPos - minPos;
r = max(sizes)/2+r_add; % radius of circle
x_centre = minPos(1) + max(sizes)/2+x_add;
y_centre = minPos(2) + max(sizes)/2+y_add;

% Data points on a perfect circle
theta = 0:pi/50:2*pi;
x_val = r * cos(theta) + x_centre;
y_val = r * sin(theta) + y_centre;

% Make two (or more) copies of this circle
moveDist = 2*r+move_add;
x_val1 = x_val - moveDist;
x_val3 = x_val + moveDist;
x_val4 = x_val - 2*moveDist;
x_val5 = x_val + 2*moveDist;
x_val6 = x_val - 3*moveDist;

y_val1 = y_val;
y_val3 = y_val;

% figure(2)
% hold on
% plot(x_val,y_val,'k','LineWidth',4)
% plot(x_val1,y_val,'k','LineWidth',4)
% plot(x_val3,y_val,'k','LineWidth',4)
% plot(x_val4,y_val,'k','LineWidth',4)
% plot(x_val5,y_val,'k','LineWidth',4)
% set(gca,'visible','off')
% axis equal

%% Plot trajectrory
clf;

% Liten: 101, 77, 81
% Medium: 4, 17, 36
% Stor: 17, 12, 35
% file1 = file;

if strcmp(exp,'large')
    file2    = 'XMLfiles/Circle_large_1agent/Stor (12)_Tracks.xml'; 
    file3    = 'XMLfiles/Circle_large_1agent/Stor (35)_Tracks.xml'; 
    kir      = abs([-0.20704; 0.49383; -3.3031]); % Stor
    img_add = 20;
elseif strcmp(exp,'medium')
    file2    = 'XMLfiles/Circle_medium_1agent/Medium (17)_Tracks.xml'; 
    file3    = 'XMLfiles/Circle_medium_1agent/Medium (36)_Tracks.xml'; 
    kir      = abs([0.18783; -0.99264; -4.6946]); % Medium
    img_add = 20;
elseif strcmp(exp,'small')
    file2    = 'XMLfiles/Circle_small_1agent/Liten (77)_Tracks.xml'; 
    file3    = 'XMLfiles/Circle_small_1agent/Liten (81)_Tracks.xml'; 
    kir      = abs([-0.2172; -1.3443; -6.272]); % Liten
    img_add = 10;
elseif strcmp(exp,'hm')

end

% file1    = file;
% file2    = [path '(12)_Tracks.xml']; 
% file3    = [path '(35)_Tracks.xml']; 
% kir      = abs([-0.20704; 0.49383; -3.3031]) % Stor


dT       = 1/25; % Time step
agent    = 1; % which agent/s to plot
digit    = 2;

[M1] = cut(file1,agent); % Turns file into a position matrix, without NaN:s
[M2] = cut(file2,agent); % Turns file into a position matrix, without NaN:s
[M3] = cut(file3,agent); % Turns file into a position matrix, without NaN:s

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

x2 = M2(agent,1,100:end);
y2 = M2(agent,2,100:end);
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
     
img = imread('hmMedBack.PNG');
img_size = r+img_add;
lightGrey1   = [0.85 0.85 0.85];

figure(1)
hold on
image('CData',img,'XData',[x_centre-img_size-moveDist x_centre+img_size-moveDist],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size x_centre+img_size],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size+moveDist x_centre+img_size+moveDist],'YData',[y_centre-img_size y_centre+img_size])

set(gca,'visible','off')
axis equal
xlim([x_centre-moveDist-img_size x_centre+moveDist+img_size])
ylim([y_centre-1.5*r y_centre+1.5*r])

p1 = plot(x1,y1,'color',c1,'DisplayName',txt1,'LineWidth',3)
p2 = plot(x2,y2,'color',c2,'DisplayName',txt2,'LineWidth',3)
p3 = plot(x3,y3,'color',c3,'DisplayName',txt3,'LineWidth',3)

plot(x_val1,y_val,'Color',lightGrey1,'LineWidth',4)
plot(x_val,y_val,'Color',lightGrey1,'LineWidth',4)
plot(x_val3,y_val,'Color',lightGrey1,'LineWidth',4)
% plot(x_val,y_val1,'k','LineWidth',6)
% plot(x_val,y_val,'k','LineWidth',4)
% plot(x_val,y_val3,'k','LineWidth',6)

set(gca,'visible','off')
axis equal
xlim([x_centre-img_size-moveDist x_centre+img_size+moveDist])
ylim([y_centre-img_size y_centre+img_size])
% xlim(x_centre-img_size x_centre+img_size]);
% ylim(y_centre-img_size y_centre+img_size]);
%hlegend = legend([p1,p2,p3],'Location','north','Interpreter','latex','FontSize',20)
%hlegend.NumColumns = 3;
% a = axes('position',get(gca,'position'),'visible','off')
% b = copyobj(a,gcf);
% c = copyobj(a,gcf);
% hlegend1 = legend(a,p1,'Location','northwest','Interpreter','latex','FontSize',20)
% hlegend2 = legend(b,p2,'Location','north','Interpreter','latex','FontSize',20)
% hlegend3 = legend(c,p3,'Location','northeast','Interpreter','latex','FontSize',20)
%% Two agents circular environment
clc;
img_add     = 20;
dT          = 1/25; % Time step
agent       = 1:2; % which agent/s to plot
digit       = 2;
values      = load('results/Lab/circle_medium_2agent.txt');
resultsPath = 'results/Lab/circle_medium_2agentSourceFiles.txt';
c_col       = 4;
chirVal     = -1;
thresh      = 0.2;
abs         = true;
val1        = -4.6;
val2        = 2;
val3        = -4.6;
fontSize    = 10;

% Agent with low chirality, most effective chirality and highest chirality,
% change abs if considering absolute values or not
% [file1,file2,file5,kir1,kir2,kir5] = findMinEffMaxGivenChir(values,resultsPath,chirVal,thresh,abs);

% Find three files where one agent has chirVal, the others are close to
% val1, val2, val3 respecitvely
[file5,file6,file3, kir5,kir6,kir3] = findSpecChirGivenChir(values,resultsPath,chirVal,thresh,val1,val2,val3);

% Min times
[file1,kir1] = findMin2Agents(values,c_col,resultsPath); % Shortest
% [file2,kir2] = findSecMin2Agents(values,c_col,resultsPath,2); % Second shortest
%[file3,kir3] = findSecMin2Agents(values,c_col,resultsPath,3);% Third shortest

% Max times
% [file1,kir1] = findMax2Agents(values,c_col,resultsPath); % Longest
[file2,kir2] = findSecMax2Agents(values,c_col,resultsPath,11); % Second longest
[file3,kir3] = findSecMax2Agents(values,c_col,resultsPath,2); % Second longest
[file4,kir4] = findSecMax2Agents(values,c_col,resultsPath,3); % Second longest
% [file3,kir3] = findSecMax2Agents(values,c_col,resultsPath,3); % Third longest

% Agent with value closest to maximum in medium circle
% [file2,kir2] = findClosestChir(values,resultsPath,chirVal); % Closest to chirVal

[x1,y1,c1,txt1,txt2,kir1] = getValues(file1,kir1);% Min time
[x2,y2,c2,txt3,txt4,kir2] = getValues(file2,kir2);% Max time
[x3,y3,c3,txt5,txt6,kir3] = getValues(file3,kir3);% Max time
[x4,y4,c4,txt7,txt8,kir4] = getValues(file4,kir4);% Max time
[x5,y5,c5,txt9,txt10,kir5] = getValues(file5,kir5);% Max time
[x6,y6,c6,txt11,txt12,kir6] = getValues(file6,kir6);% Max time


% Position by moving with moveDist
x1 = x1-moveDist;
%x2 = x2;
x3 = x3+moveDist;
x4 = x4+2*moveDist;
x5 = x5-2*moveDist;
x6 = x6-3*moveDist;

% Load background picture
img = imread('hmMedBack.PNG');
img_size = r+img_add;
lightGrey1   = [0.85 0.85 0.85];

figure(1)
hold on
% Plots background image
image('CData',img,'XData',[x_centre-img_size-moveDist x_centre+img_size-moveDist],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size x_centre+img_size],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size+moveDist x_centre+img_size+moveDist],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size-2*moveDist x_centre+img_size-2*moveDist],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size+2*moveDist x_centre+img_size+2*moveDist],'YData',[y_centre-img_size y_centre+img_size])
image('CData',img,'XData',[x_centre-img_size-3*moveDist x_centre+img_size-3*moveDist],'YData',[y_centre-img_size y_centre+img_size])

% Plot trajectories
p1 = plot(x1(:,1),y1(:,1),'color',c1(1,:),'DisplayName',txt1,'LineWidth',3);
p2 = plot(x1(:,2),y1(:,2),'color',c1(2,:),'DisplayName',txt2,'LineWidth',3);
p3 = plot(x2(:,1),y2(:,1),'color',c2(1,:),'DisplayName',txt3,'LineWidth',3);
p4 = plot(x2(:,2),y2(:,2),'color',c2(2,:),'DisplayName',txt4,'LineWidth',3);
p5 = plot(x3(:,1),y3(:,1),'color',c3(1,:),'DisplayName',txt5,'LineWidth',3);
p6 = plot(x3(:,2),y3(:,2),'color',c3(2,:),'DisplayName',txt6,'LineWidth',3);

p7 = plot(x4(:,1),y4(:,1),'color',c4(1,:),'DisplayName',txt7,'LineWidth',3);
p8 = plot(x4(:,2),y4(:,2),'color',c4(2,:),'DisplayName',txt8,'LineWidth',3);
p9 = plot(x5(:,1),y5(:,1),'color',c5(1,:),'DisplayName',txt9,'LineWidth',3);
p10 = plot(x5(:,2),y5(:,2),'color',c5(2,:),'DisplayName',txt10,'LineWidth',3);
p11 = plot(x6(:,1),y6(:,1),'color',c6(1,:),'DisplayName',txt11,'LineWidth',3);
p12 = plot(x6(:,2),y6(:,2),'color',c6(2,:),'DisplayName',txt12,'LineWidth',3);

% Plot collisions and find times
[time1] = checkCollision(x1,y1);
[time2] = checkCollision(x2,y2);
[time3] = checkCollision(x3,y3);
[time4] = checkCollision(x4,y4);
[time5] = checkCollision(x5,y5);
[time6] = checkCollision(x6,y6);

% Fake points used for time legend
f1 = plot(10000,10000,'w','DisplayName',time1)
f2 = plot(10000,10000,'w','DisplayName',time2)
f3 = plot(10000,10000,'w','DisplayName',time3)
f4 = plot(10000,10000,'w','DisplayName',time4)
f5 = plot(10000,10000,'w','DisplayName',time5)
f6 = plot(10000,10000,'w','DisplayName',time6)

% Plots circles
plot(x_val1,y_val,'Color',lightGrey1,'DisplayName',time1,'LineWidth',4);
plot(x_val,y_val,'Color',lightGrey1,'DisplayName',time2,'LineWidth',4);
plot(x_val3,y_val3,'Color',lightGrey1,'DisplayName',time3,'LineWidth',4);
plot(x_val4,y_val,'Color',lightGrey1,'DisplayName',time4,'LineWidth',4);
plot(x_val5,y_val3,'Color',lightGrey1,'DisplayName',time5,'LineWidth',4);
plot(x_val6,y_val3,'Color',lightGrey1,'DisplayName',time6,'LineWidth',4);

% Axis settings and legends
set(gca,'visible','off');
axis equal;
xlim([x_centre-3*moveDist-img_size x_centre+3*moveDist+img_size]);
ylim([y_centre-img_size y_centre+img_size]);
%hlegend = legend([p1,p2,p3],'Location','north','Interpreter','latex','FontSize',20)
%hlegend.NumColumns = 3;
a = axes('position',get(gca,'position'),'visible','off');
b = copyobj(a,gcf);
c = copyobj(b,gcf);
g = copyobj(c,gcf);
h = copyobj(c,gcf);
i = copyobj(b,gcf);
dist = 0.11
% hlegend1 = legend(a,[p1,p2],'Location','north','Interpreter','latex','FontSize',fontSize);
% hlegend2 = legend(b,[p3,p4],'Location','north','Interpreter','latex','FontSize',fontSize);
% hlegend3 = legend(c,[p5,p6],'Location','north','Interpreter','latex','FontSize',fontSize); 
% hlegend7 = legend(g,[p7,p8],'Location','north','Interpreter','latex','FontSize',fontSize);
% hlegend8 = legend(h,[p9,p10],'Location','north','Interpreter','latex','FontSize',fontSize);
% hlegend9 = legend(i,[p11,p12],'Location','north','Interpreter','latex','FontSize',fontSize);
hl = legend(i,[p11,p12,p9,p10,p1,p2,p3,p4,p5,p6,p7,p8],'Location','north','Interpreter','latex','FontSize',fontSize)
hl.NumColumns = 6
annotation('textbox',[0.14 0.3 0.1 0.1],'String',time6,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.14+dist 0.3 0.1 0.1],'String',time5,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.14+2*dist+0.002 0.3 0.1 0.1],'String',time1,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.14+3*dist+0.007 0.3 0.1 0.1],'String',time2,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.14+4*dist+0.009 0.3 0.1 0.1],'String',time3,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
annotation('textbox',[0.14+5*dist+0.01 0.3 0.1 0.1],'String',time4,'FitBoxToText','on','Interpreter','latex','FontSize',fontSize);
hold off

%% Homogenous, find files
%clear; clc; clf;
%Image dimensions
% Film 43: min y
% Film 44: min x
% FIlm 74: max x
% Film 72: max y
agent   = 1;
c       = load('results/Lab/hm1agent.txt');
digit   = 2;
exp_path = 'XMLfiles/Homogen_1agent/';
results_path = 'results/Lab/hm1agentSourceFiles.txt';
sf = 1;

% minx, max_X
f1 = [exp_path '044_Tracks.xml'];
f2 = [exp_path '074_Tracks.xml'];
[P1] = cut(f1,agent); % Turns file into a position matrix, without NaN:s
[P2] = cut(f2,agent); % Turns file into a position matrix, without NaN:s
% min_y, max_y
f3 = [exp_path '043_Tracks.xml'];
f4 = [exp_path '072_Tracks.xml'];
[P3] = cut(f3,agent); % Turns file into a position matrix, without NaN:s
[P4] = cut(f4,agent); % Turns file into a position matrix, without NaN:s

x_min_pos = min(P1(:,1,:))*sf;
x_max_pos = max(P2(:,1,:))*sf;
x_centre = (x_max_pos+x_min_pos)/2;

y_min_pos = min(P3(:,2,:))*sf;
y_max_pos = max(P4(:,2,:))*sf;
y_centre = (y_max_pos+y_min_pos)/2;

img = imread('medBack.PNG');
img_size = [(x_max_pos-x_min_pos)/2 (y_max_pos-y_min_pos)/2];

moveDist = 2*img_size(1)+100;

kir_abs = abs(c(:,1));
[min_kir, ind_min] = min(kir_abs)
[max_kir, ind_max] = max(kir_abs)

area = c(:,2);
comp_val = ones(length(area),1);
difference = abs(comp_val-area);
[mid_area_diff, ind_mid] = min(difference);
mid_area = area(ind_mid);
mid_kir = kir_abs(ind_mid);

% min_kir = kir_abs(3);
% mid_kir = kir_abs(11);
min_kir = kir_abs(47);
kir = [min_kir mid_kir max_kir];

% Find files
file1 = 'XMLfiles/Homogen_1agent/072_Tracks.xml'; %index 47

% fileID = fopen(results_path,'r');
% for i = 1:ind_min
%     file1 = fgetl(fileID);
% end
% fclose(fileID);
% file1 = 'XMLfiles/KorsTefat_1agent/0.54_Tracks.xml'
% fileID = fopen(results_path,'r');
% for i = 1:ind_mid
%     file2 = fgetl(fileID);
% end
% fclose(fileID);
%file2 = 'XMLfiles/KorsTefat_1agent/0.76_Tracks.xml'
file2 = 'XMLfiles/HomogenLeft_1agent/6_Tracks.xml';
fileID = fopen(results_path,'r');
for i = 1:ind_max
    file3 = fgetl(fileID);
end

[M1,length1,times1] = cut(file1,agent); % Turns file into a position matrix, without NaN:s
[M2,length2,times2] = cut(file2,agent); % Turns file into a position matrix, without NaN:s
[M3,length3,times3] = cut(file3,agent); % Turns file into a position matrix, without NaN:s

c1 = chir2color(abs(kir(1)));
c2 = chir2color(abs(kir(2)));
c3 = chir2color(abs(kir(3)));

x1 = M1(agent,1,1:end-1)*sf;
y1 = M1(agent,2,1:end-1)*sf;
x1 = x1(:,:)'-moveDist;
y1 = y1(:,:)';

x2 = M2(agent,1,:)*sf;
y2 = M2(agent,2,:)*sf;
x2 = x2(:,:)';
y2 = y2(:,:)';

x3 = M3(agent,1,:)*sf;
y3 = M3(agent,2,:)*sf;
x3 = x3(:,:)'+moveDist;
y3 = y3(:,:)';

num = round(kir,digit,'significant');


txt1 = ['$\omega$ = ',num2str(num(1)),' rad/s']; 
txt2 = ['$\omega$ = ',num2str(num(2)),' rad/s'];
txt3 = ['$\omega$ = ',num2str(num(3)),' rad/s'];

figure(1)
hold on
image('CData',img,'XData',[x_centre-img_size(1) x_centre+img_size(1)],'YData',[y_centre-img_size(2) y_centre+img_size(2)])
image('CData',img,'XData',[x_centre-img_size(1)-moveDist x_centre+img_size(1)-moveDist],'YData',[y_centre-img_size(2) y_centre+img_size(2)])
image('CData',img,'XData',[x_centre-img_size(1)+moveDist x_centre+img_size(1)+moveDist],'YData',[y_centre-img_size(2) y_centre+img_size(2)])

p1 = plot(x1,y1,'color',c1,'DisplayName',txt1,'LineWidth',3);
p2 = plot(x2,y2,'color',c2,'DisplayName',txt2,'LineWidth',3);
p3 = plot(x3,y3,'color',c3,'DisplayName',txt3,'LineWidth',3);

set(gca,'visible','off')
axis equal
xlim([x_centre-moveDist-img_size(1) x_centre+moveDist+img_size(1)])
ylim([y_centre-img_size(2) y_centre+img_size(2)])

a = axes('position',get(gca,'position'),'visible','off');
b = copyobj(a,gcf);
c = copyobj(a,gcf);
legend(a,p1,'Location','northwest','Interpreter','latex','FontSize',20)
legend(b,p2,'Location','north','Interpreter','latex','FontSize',20)
legend(c,p3,'Location','northeast','Interpreter','latex','FontSize',20)

hold off

%% Find file
function file = findFile(ind,path)
    fileID = fopen(path,'r');
    for i = 1:ind
        file = fgetl(fileID); % xml-file for 
    end
    fclose(fileID);
end

%% Do stuff with file

% get x & y
% color
% 

function [x,y,col,txt1,txt2,orderedKir] = getValues(file,kir)
digit = 2;
agent = 1:2

[M1] = cut(file,agent); % Turns file into a position matrix, without NaN:s

col = chir2color(abs(kir));

x = M1(agent,1,:);
y = M1(agent,2,:);
x = x(:,:)';
y = y(:,:)';

%Check order of agents
if y(1,1)<y(1,2)
   tempx = x(:,1);
   x(:,1) = x(:,2);
   x(:,2) = tempx;
   tempy = y(:,1);
   y(:,1) = y(:,2);
   y(:,2) = tempy;
   tempcol = col(1,:)
   col(1,:) = col(2,:)
   col(2,:) = tempcol;
   tempKir = kir(1);
   kir(1) = kir(2);
   kir(2) = tempKir;
end
orderedKir = kir;

num = round(orderedKir,digit,'significant');

txt1 = ['$\omega$ = ',num2str(num(1)),' rad/s$\:$'];
txt2 = ['$\omega$ = ',num2str(num(2)),' rad/s$\:$']; 


end

%% Find max time_file
function [file,kir] = findMax2Agents(c,c_col,resultsPath)

% Which c value to consider
times = c(:,c_col);
% Finds maximum value in vector and its index
[t_max, ind_max] = max(times);

kir = [c(ind_max,1) c(ind_max+1,1)];

file = findFile(ceil(ind_max/2),resultsPath);
end
%% Find min time file
function [file,kir] = findMin2Agents(c,c_col,resultsPath)

% Which c value to consider
times = c(:,c_col);
% Finds maximum value in vector and its index
[t_min, ind_min] = min(times);

kir = [c(ind_min,1) c(ind_min+1,1)];

file = findFile(ceil(ind_min/2),resultsPath);
end

%% Find order of min file
function [file,kir] = findSecMin2Agents(c,c_col,resultsPath,num)

num = num*2-1;
% Which c value to consider
times = c(:,c_col);
[ordered_times, ind] = sort(times);
% Finds maximum value in vector and its index
t_val = ordered_times(num);

orig_ind = ind(num);
kir = [c(orig_ind,1) c(orig_ind+1,1)];

file = findFile(ceil(orig_ind/2),resultsPath);

end

%% Find order of max file
function [file,kir] = findSecMax2Agents(c,c_col,resultsPath,num)

num = length(c(:,c_col))-num*2+1;
% Which c value to consider
times = c(:,c_col);
[ordered_times, ind] = sort(times);
% Finds maximum value in vector and its index

t_val = ordered_times(num);

orig_ind = ind(num);
kir = [c(orig_ind,1) c(orig_ind+1,1)];

file = findFile(ceil(orig_ind/2),resultsPath);

end
%% Find Agent w chirality close to something

function [file,kir,ind_minDiff] = findClosestChir(c,resultsPath,chirVal)
    
    chiralities = c(:,1);
    comp_val = ones(length(chiralities),1)*chirVal;
    difference = abs(comp_val-chiralities);
    [chir_diff, ind_minDiff] = min(difference);
    closestChir = chiralities(ind_minDiff);
    
    if (isItEven(ind_minDiff))
        % Even number, i.e second in a pair
        ind_minDiff = ind_minDiff-1;
    end
    
    kir = [c(ind_minDiff,1) c(ind_minDiff+1,1)];
    file = findFile(ceil(ind_minDiff/2),resultsPath);
    
end
%% Is it even?
function [even] = isItEven(num)

    if (num/2-floor(num/2) == 0)
        % Even number, i.e second in a pair
        even = 1;
    else
        even = 0;
    end
end

%% Find min kir, max kir and most effective kir for a certain chirality chirVal

function [file1,file2,file3,kir1,kir2,kir3] = findMinEffMaxGivenChir(c,resultsPath,chirVal,thresh,abs_val)

kirs = c(:,1);
times = c(:,4);

ind = zeros(1,length(kirs));
pos = 1;
for i=1:length(kirs)
   if (abs(kirs(i)-chirVal)<thresh)
       if (isItEven(i))
           save = i-1
       else
           save = i;
       end
       ind(pos) = save;
       pos=pos+1;
   end
end
nonzeroindices = find(nonzeros(ind))
ind = ind(nonzeroindices)

relevant_kirs = kirs(ind)
relevant_times = times(ind)

% Find a min kir, max kir and most effective kir
if (abs_val == true)
    [kir_min, ind_min] = min(abs(relevant_kirs));
    [kir_max, ind_max] = max(abs(relevant_kirs));
    kir_min = relevant_kirs(ind_min);
    kir_max = relevant_kirs(ind_max);
else
    [kir_min, ind_min] = min(relevant_kirs);
    [kir_max, ind_max] = max(relevant_kirs);
 end
% Most effective
[time_min, ind_short] = min(relevant_times);
kir_eff = relevant_kirs(ind_short);

% Find original indices in saved data
indices = [ind(ind_min) ind(ind_short) ind(ind_max)];
other_kir1 = 1;
other_kir2 = 1;
other_kir3 = 1;
if isItEven(indices(1))
    indices(1) = indices(1)-1
    other_kir1 = 0;
    %kir1 = [kirs(indices(1)) kirs(indices(1)+1)]
end
if isItEven(indices(2))
    indices(2) = indices(2)-1
    other_kir2 = 0;
    %kir2 = [kirs(indices(2)) kirs(indices(2)+1)]
end
if isItEven(indices(3))
    indices(3) = indices(3)-1
    other_kir3 = 0;
    %kir3 = [kirs(indices(3)) kirs(indices(3)+1)]
end

%kirs = [kir_min kir_eff kir_max];
kir1 = [kir_min kirs(indices(1)+other_kir1)];
kir2 = [kir_eff kirs(indices(2)+other_kir2)];
kir3 = [kir_max kirs(indices(3)+other_kir3)];
file1 = findFile(ceil(indices(1)/2), resultsPath);
file2 = findFile(ceil(indices(2)/2), resultsPath);
file3 = findFile(ceil(indices(3)/2), resultsPath);

end

%% Find three specific chiralities for one agent having given chirality chirVal
function [file1,file2,file3,kir1,kir2,kir3] = findSpecChirGivenChir(c,resultsPath,chirVal,thresh,val1,val2,val3)

kirs = c(:,1);
times = c(:,4);

ind = zeros(1,length(kirs));
pos = 1;
for i=1:length(kirs)
   if (abs(kirs(i)-chirVal)<thresh)
       if (isItEven(i))
           save = i-1
       else
           save = i;
       end
       ind(pos) = save;
       pos=pos+1;
   end
end
nonzeroindices = find(nonzeros(ind))
ind = ind(nonzeroindices)

relevant_kirs = kirs(ind)
relevant_times = times(ind)

% Close to val1 val2 val3
[f1,kir1,ind1] = findClosestChir(relevant_kirs,resultsPath,val1)
[f2,kir2,ind2] = findClosestChir(relevant_kirs,resultsPath,val2)
[f3,kir3,ind3] = findClosestChir(relevant_kirs,resultsPath,val3)

% Find original indices in saved data
indices = [ind(ind1) ind(ind2) ind(ind3)];
other_kir1 = 1;
other_kir2 = 1;
other_kir3 = 1;
if isItEven(indices(1))
    indices(1) = indices(1)-1
    other_kir1 = 0;
    %kir1 = [kirs(indices(1)) kirs(indices(1)+1)]
end
if isItEven(indices(2))
    indices(2) = indices(2)-1
    other_kir2 = 0;
    %kir2 = [kirs(indices(2)) kirs(indices(2)+1)]
end
if isItEven(indices(3))
    indices(3) = indices(3)-1
    other_kir3 = 0;
    %kir3 = [kirs(indices(3)) kirs(indices(3)+1)]
end

%kirs = [kir_min kir_eff kir_max];
kir1 = [kir1(1) kirs(indices(1)+other_kir1)];
kir2 = [kir2(1) kirs(indices(2)+other_kir2)];
kir3 = [kir3(1) kirs(indices(3)+other_kir3)];
file1 = findFile(ceil(indices(1)/2), resultsPath);
file2 = findFile(ceil(indices(2)/2), resultsPath);
file3 = findFile(ceil(indices(3)/2), resultsPath);


end

%%
function [time] = checkCollision(x,y)
dT = 0.04;
if (length(x)<1380)
    midP = [(x(end,1)+x(end,2))/2, (y(end,1)+y(end,2))/2];
    plot(midP(1),midP(2),'kx','LineWidth',3,'MarkerSize',20);
    time = [num2str(length(x)*dT) ' sekunder$\quad$'];
else
    time = 'Maximal tid$\quad$';
end
end