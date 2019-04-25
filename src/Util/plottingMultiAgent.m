%plottingMultiAgent


%% Create plots
load('results/Final_results/2019411-1611_multiAgentCircle_R1_N100_w75.mat')

X = w;
Y = w;
Z_1 = meanTotalTime + meanTotalTime' - diag(diag(meanTotalTime));
Z_2 = 1./Z_1;


figure(1)
surf(X,Y,Z_1)
title("Time")

figure(2)
surf(X,Y,Z_2)
title("Efficiency")

Z_3 = percentMax + percentMax' - diag(diag(percentMax));

figure(3)
surf(X,Y,Z_3)
title("Fraction failures")

Z_4 = Z_3;

for i = 1:size(Z_4,1)
    for j = 1:size(Z_4,2)
        if Z_4(i,j) ~= 1 && Z_4(i,j) ~= 0
            Z_4(i,j) = 0.5;
        end
    end
end

figure(4)
surf(X,Y,Z_4)
title("Win and fail")

%%Formatting 1

ax_Font = 40;           % Fontsizes on axis
tit_Font = 30;          % Fontsizes on titles

figure(1)

view(2)

axis('square')
title('Tid till m�te f�r tv� agenter', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Kiralitet agent A (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet agent B (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
bar = colorbar;
set(get(bar,'label'),'string','Tid (s)', 'Interpreter', 'latex', 'fontsize', ax_Font); 

%%Formatting 2

figure(2)

view(2)

axis('square')
title('Invers av tid till m�te av tv� agenter (Effektivitet)', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Kiralitet agent A (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet agent B (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
bar = colorbar;
set(get(bar,'label'),'string','Effektivitet (s$^{-1}$)', 'Interpreter', 'latex', 'fontsize', ax_Font);

%%Formatting 3

figure(3)

view(2)

axis('square')
title('Andel agentpar som inte hittar varandra vid en maxtid p� 300 sekunder', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Kiralitet agent A (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet agent B (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
bar = colorbar;
set(get(bar,'label'),'string','Andel', 'Interpreter', 'latex', 'fontsize', ax_Font);
