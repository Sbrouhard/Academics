% Chris Danforth, 1/10/23
% example of finding a periodic orbit of the logistic map g(x) = ax(1-x)

clear all		% always start with these 2 lines
close all

format long		% print all sig figs unless instructed otherwise


%Note: I found these orbits by placing the entire function into a loop with
% a value a_iterate. This is what I'm using to slightly push the initial
% condition
% higher. The value 3.569945261497732 which is the leading term in the
% "pushing forward" was found by starting at ~3.6 and just continuing to
% find higher orbits and replacing the current leading term with the new "working" one. 
% As I started running
% too far forward and diverged, I added another 0 to the multiplier of the a_iterate and
% re-ran the program, to make the increases smaller. I eventually landed on
% this giving a 2048 orbit, and called it there!


a_max_value = 0;
a_max_orbit = 0;

for a_iterate = 1:1000
    %% initial stuff
    itotal = 10^8;		% total number of iterates
    x = zeros(itotal,1);	% fill orbit with zeros
    x(1) = rand;		% random initial condition		% parameter of logistic map, change this to get different sinks, try other values listed
    a =  (3.569945261497732 + sqrt((0.0000000000000001 * a_iterate)));	% parameter of logistic map, change this to get different sinks, try other values listed
    
    %% iterate the map for a long time
    for i = 1:itotal-1
      x(i+1) = a*x(i)*(1-x(i));		% logistic map, g(x)=ax(1-x)
    end
    
    T = 35;
    %% print the last T points of the orbit
    disp(' ')
    disp('  Iterate        X')
    disp('---------------------------------------------')
    for i = itotal-T+1:itotal					% only showing T points on the orbit
      disp(sprintf('   %10.0f  %1.14f',i,x(i)))
    end
    disp(' ')
    
    %% determine if there is a periodic orbit of period less than 2^S
    S = 27;					% length of longest periodic orbit (log_2) we can find, maybe?
    tol = 10^(-8);				% tolerance within which I believe I've found a repeat
    inside=0;				% trigger telling me if I've been inside the 'if' loop below
    for i = 0:S				% loop looking for orbit of period 2^(i)
      difference(i+1) = norm(x(end)-x(end-2^i));		% compute difference between orbit points
      if (difference(i+1) < tol) & (inside == 0)		% am I within tolerance? have I already qualified?
        disp(sprintf('This orbit repeats every %s iterates',int2str(2^i)))
        if (2^i > a_max_orbit)
            a_max_orbit = 2^i;
            a_max_value = a
        end
        inside=1;		% don't come back inside this loop, see what happens without this
      end
    end
   disp(sprintf('The a value is %1.14f', a))
end
disp(sprintf('Final Max Orbit: %s   Final max value for a: %1.16f', int2str(a_max_orbit), a_max_value))


