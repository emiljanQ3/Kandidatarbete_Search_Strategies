%% Find circle

file = 'XMLfiles/Circle_small_1agent/Liten (101)_Tracks.xml'; 
agent = 1; % which agent/s to plot
[pos_a,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

maxPos = max(max(pos_a,[],1),[],3);
minPos = min(min(pos_a,[],1),[],3);

size = maxPos - minPos;
r = max(size)/2 + 2;
x_centre = minPos(1) + max(size)/2;
y_centre = minPos(2) + max(size)/2;

theta = 0:pi/50:2*pi;
x_val = r * cos(theta) + x_centre;
y_val = r * sin(theta) + y_centre;


%% Plot trajectrory
clf;
expName = 'circle_small_1agent';
name = join(['results/Lab/' expName '.txt']);
c = load(name);
all_kir = c(:,1);

dT = 1/25; % Time step
[all_kir_sorted, sortOrder] = sort(all_kir);
c1 = chir2color(all_kir_sorted);

file = 'XMLfiles/Circle_small_1agent/Liten (65)_Tracks.xml'; 
agent = 1; % which agent/s to plot
[pos_a,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

myCircle = [pos_a(1,1,1:1), pos_a(1,2,1:1)];
[r, indice] = splitPositionDataPartitioned(pos_a,700, myCircle);
[kir,D_r,v] = getKompSpiral(r,dT,1,6,60);

i = 1;
while kir ~= all_kir_sorted(i)
   i = i+1; 
end

% SourceFiles = join(['results/Lab/' expName 'SourceFiles.txt']);
% fileID = fopen(name,'r');
% A = fread(fileID,'*uint8');
% A_text = sprintf('%.2x',A);

x = pos_a(agent,1,:);
y = pos_a(agent,2,:);
x = x(:,:)';
y = y(:,:)';

figure(1)
hold on
plot(x_val,y_val,'k','LineWidth',5)
plot(x,y,'color',all_kir_sorted(i))
axis equal