% Makes a gif of a simulation
clear;clc;clf;
load('results/SimulatedTracks/Homogen/infinite_w0.1557.mat')

path        = ['Gifs/Homogen', num2str(w), '/'];
status      = mkdir(path);
headSize    = 17;
extraSize   = 3;

maxPos = max(max(pos_a,[],1),[],3)+extraSize*ones(1,2);
minPos = min(min(pos_a,[],1),[],3)-extraSize*ones(1,2);

figure(1)
hold on
set(gca,'visible','off')
xlim([minPos(1) maxPos(1)]);
ylim([minPos(2) maxPos(2)]);
%axis equal

h2 = plot(pos_a(:,1,1),pos_a(:,2,1),'Marker','.','MarkerSize',headSize,'color','k')

for time = 2:numTimeSteps
    pause(dT)
    delete(h2)
    h =  plot([pos_a(:,1,time),pos_a(:,1,time-1)]',[pos_a(:,2,time),...
         pos_a(:,2,time-1)]','color',[0.4660 0.6740 0.1880],...
         'LineWidth',3);
    h2 = plot(pos_a(:,1,time),pos_a(:,2,time),'Marker','.','MarkerSize',headSize,'color','k')
    saveas(gcf,[path, num2str(time), '.png']);
end
