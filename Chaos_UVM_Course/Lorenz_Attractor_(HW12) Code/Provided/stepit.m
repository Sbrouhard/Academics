% stepit.m
%
% For use with a system of autonomous differential equations i.e. dy/dt=f(y)
%
% Performs numerical integration of a function by Runge-Kutta
% method, over the range (t0, t0 + step). Therefore
% only one iteration is performed.
%
% Useage;
% stepit('func', x0, parameters, stepsize)
%
% where func is the name of the function to be integrated, funcarg is the 
% left hand limit of the integration, parameters contains any necessary
% parameters to pass to the function and stepsize is the 
% integration timestep.
%
% You would normally include this in a loop.
%
% Modified by Nick Allgaier 4/27/ll

function x = stepit(func, x, p, step)
	s1 = feval(func,x,p);
	s2 = feval(func,x + step*s1/2,p);
	s3 = feval(func,x + step*s2/2,p);
	s4 = feval(func,x + step*s3,p);
	x = x + step*(s1 + 2*s2 + 2*s3 + s4)/6;