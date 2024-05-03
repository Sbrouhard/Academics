%% Plot the Henon Map
%% Chris Danforth 2/1/23

clear all
close all
format long		% show me lots of digits

%% initialize parameters
%a = 1.28; b = -0.3;	% smooth boundary fig 2.3
%a = 1.4; b = -0.3;	% fractal boundary fig 2.3
%a = 2; b = -0.3;	% attractor fig 2.11
a = 1.950300; b = -0.3;		% smooth boundary fig 2.9, 2.10
%a = 0.43; b = 0.4;

%% initialize plot
figure
xmin = -2.5; xmax  = 2.5;
ymin = -2.5; ymax  = 2.5;
colormap([0 0 0;1 1 1])
axis([xmin xmax -ymax ymax])
fsize=16;
fsizea=30;
set(gca,'xtick',[xmin:1:xmax],'FontSize',fsize)
set(gca,'ytick',[-ymax:1:ymax],'FontSize',fsize)
xlabel('x','FontSize',fsizea); ylabel('y','FontSize',fsizea)
title('Henon Map','FontSize',fsizea)
rect = [0,30,560*1.0,560*0.8];
set(gcf,'Position',rect);
set(gcf,'PaperPositionMode','auto');   

%% ask for an initial condition
disp('Click on a point and press enter');
[x,y] = ginput;
disp(sprintf('You chose x = %0.3f, y = %0.3f',x,y)); pause(1)
hold on
iterates = 2000;		% how many times we iterate the map
lag = 50;		% how many iterates are plotted on figure at any time

for j = 1:iterates
  if mod(j,10) == 0
    disp(sprintf('Iterate %s',int2str(j)))
  end
  if j == lag
    disp('Removing Old Iterates')
  end
  yold = y; y=a1x + b1y; x=a-x^2+b*yold;	% this is the henon map
  disp(sprintf('(x,y) coordinates: %d,%d',x,y))
  h(j)=plot(x,y,'r+','MarkerSize',15);	% plot the current orbit point
  if j > lag
    delete(h(j-lag))			% remove old iterates
  end
  pause(.1)				% pause for drama, remove for insane speed
end
