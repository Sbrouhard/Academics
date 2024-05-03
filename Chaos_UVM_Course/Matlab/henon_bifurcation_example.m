%% Make bifurcation diagram for Henon Map
%% Chris Danforth 2/17/2022
import java.util.*
clear all
close all

%% define parameters
%amin = 1.925; amax = 1.975; b = -0.3; 
amin = 0; amax = 2; b = -0.3; 
N = 10^6;
%N = 10^4;
da = 10^-3; IC = [0;2];
x = zeros(N,1); y = zeros(N,1);

%% make empty figure 
figure; 

%% iterate over parameter values a
%for a = linspace(amin,amax,25)

previous_period = 1;
alphavalues = Stack();

for a = amin:da:amax
  disp(a)
  x(1) = IC(1);
  y(1) = IC(2);
  
  %% iterate the map for this a value
  for i = 1:N
    x(i+1) = a-x(i)^2+b*y(i);
    y(i+1) = x(i);
  end

  %% determine if there is a periodic orbit of period less than 2^S
  %% we do this to avoid plotting thousands of points on the same pixel (e.g. converged to a fixed point sink)
  S = 9;                                 % length of longest periodic orbit (log_2) we can find, maybe?
  p = S;				  % initialize length of the orbit we will plot
  tol = 10^(-5);                         % tolerance within which I believe I've found a repeat
  inside=0;                               % trigger telling me if I've been inside the 'if' loop below
  for i = 0:S                             % loop looking for orbit of period 2^(i)
    difference(i+1) = norm(x(end)-x(end-2^i));            % compute difference between orbit points
    if (difference(i+1) < tol) & (inside == 0)            % am I within tolerance? have I already qualified?
      inside=1;           % don't come back inside this loop, see what happens without this
      p = i;		  % save the period we discovered, we will only plot the last 2^p points
    end
  end
  disp(2^p)

  if(p > previous_period && p ~= 0)
      disp("PREV and CURRENT")
      fprintf("previous is %.6f and current is %.6f \n", previous_period, p)
      previous_period=p;
      stringtoshow = sprintf("%.6f at alpha value %.6f:) + ", p, a);
      alphavalues.push(stringtoshow);
  end

   % plot the final 2^p iterates for this a value
   plot(a,x(end-2^p:end),'k.','MarkerSize',1);
   % hold on;					% don't erase old points

  
  
end

disp(alphavalues);

 %% spruce up the plot
 % xlabel('parameter a'); ylim([-0.5 2]); xlim([amin amax]); 
 % h = get(gca, 'xlabel'); set(h, 'FontSize', 32)
 % h = get(gca, 'ylabel'); set(h, 'FontSize', 32)
 % set(gca, 'FontSize', 24); rect = [0,0,560*1.7,560*1.7];
 % set(gcf,'Position',rect); set(gcf,'PaperPositionMode','auto'); set(gca,'XMinorTick','on');% 
 % %% print it to a file I can see
 % print('-dtiff','bifurcation_1')
