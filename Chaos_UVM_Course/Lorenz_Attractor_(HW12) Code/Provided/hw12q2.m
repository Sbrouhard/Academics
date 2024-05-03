%% Nick Allgaier 4/27/11
%% Code to calculate LE for Lorenz 63

clear all 
close all

%% Lorenz system
f = @(x,par) [par(1)*(x(2)-x(1));
              par(2)*x(1)-x(2)-x(1)*x(3);
              x(1)*x(2)-par(3)*x(3)];
%% Jacobian
Df = @(x,par) [-par(1),par(1),0;
		par(2)-x(3),-1,-x(1);
		x(2),x(1),-par(3)];

%% Standard parameters
par = zeros(3,1);
s=10; par(1)=s;
r=28; par(2)=r;
b=8/3; par(3)=b;

%% IC
x0=[1.508870;-1.531271;25.46091];	%  IC 
N = 10^4;
x = zeros(3,N); x(:,1) = x0;
for i=2:N, x(:,i) = stepit(f,x(:,i-1),par,0.01); end

%% Compute Lyapunov exponents
altpar = par;
altpar(2) = 12;
LEcalc(f,Df,altpar,x(:,end),'ode') % for r=12
% altpar(2) = 24.5;
% LEcalc(f,Df,altpar,x(:,end),'ode') % for r=24.5
% LEcalc(f,Df,par,x(:,end),'ode')    % for r=28

