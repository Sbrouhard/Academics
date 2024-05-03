clear all
close all

%% initialize parameters
global b s r		
b = 8/3; s = 10; r = 28;
tspan = [0,10^3]; tstep = 0.01;
x0 = [-12.6480  -13.9758   30.9758];

%% lorenz equations
lorenz = @(t,x) [s*(x(2)-x(1));r*x(1)-x(2)-x(1)*x(3);x(1)*x(2)-b*x(3)];

%% integrate model
[t,y] = rk2(lorenz,tspan,x0,tstep);

%% make plot
plot3(y(:,1),y(:,2),y(:,3),'k.','MarkerSize',1);  grid on;
disp('Grab the arrow and rotate the attractor')

figure
boxcount3(y,6,1)
