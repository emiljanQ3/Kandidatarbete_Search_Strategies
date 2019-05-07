%Testing area elements
L = 0.156/1.7;
R = 1/1.7;

plot(R*cos(linspace(0,2*pi)),R*sin(linspace(0,2*pi)))
grid on
xticks(-R:L:R);
yticks(-R:L:R);