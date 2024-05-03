%% Plot the Mandelbrot set and Julia Set with interactive orbit generation
%% Not commented particularly well
%% Chris Danforth 2/28/07

clear all
close all

%% first the mandelbrot set
Nmax = 30; 
scale = 0.003;
%scale = 0.01;
xmin = -2.4; xmax  = 1.2;
ymin = -1.5; ymax  = 1.5;

disp('Be patient, I am working')

%% Generate x and y coordinates and z complex values
[x,y]=meshgrid(xmin:scale:xmax,ymin:scale:ymax);
z = x+1i*y;

%% Generate w accumulation matrix and k counting matrix
w = zeros(size(z));
k = zeros(size(z));

N = 0;
while N<Nmax & ~all(k(:))
    w = w.^2+z;
    N = N+1;
    k(~k & abs(w)>4) = N;
end
k(k==0) = Nmax;
figure
s = pcolor(x, y, mod(k, 2));
colormap([0 0 0;1 1 1])
set(s,'edgecolor','none')
axis([xmin xmax -ymax ymax])
fsize=16;
fsizea=30;
set(gca,'xtick',[xmin:0.4:xmax],'FontSize',fsize)
set(gca,'ytick',[-ymax:0.5:ymax],'FontSize',fsize)
xlabel('Real(z)','FontSize',fsizea)
ylabel('Imaginary(z)','FontSize',fsizea)
title('Mandelbrot Set','FontSize',fsizea)
rect = [0,30,1000,850];
%rect = [0,30,560*1.0,560*0.8];
%rect = [0,0,560*1.0,560*0.8];
set(gcf,'Position',rect);
set(gcf,'PaperPositionMode','auto');   


% now the julia set
disp('Click on a point and              THEN PRESS ENTER'); pause(1)
disp('The point you choose will define c in the map f(z) = z^2 + c'); pause(1)
disp('Hint:  A point /inside/ the set will result in more interesting stuff'); pause(1)
disp('Another Hint:  Try different /lobes/ around the main lobe')
[xup,yup] = ginput;
c = complex(0,1.1);
%c = complex(0.32,0.043);		% interesting period 11 sink
disp(sprintf('You chose c = %0.3f + %0.3fi',xup,yup)); pause(1)
disp('Be patient, I am working')

Nmax = 50; 
%scale = 0.001;
scale = 0.003;
xmin = -2; xmax  = 2;
ymin = -2; ymax  = 2;

% Generate x and y coordinates and z complex values
[x,y]=meshgrid(xmin:scale:xmax,ymin:scale:ymax);
z = x+1i*y;

% Generate w accumulation matrix and k counting matrix
w = zeros(size(z));
k = zeros(size(z));

N = 0;
while N<Nmax & ~all(k(:))
    z = z.^2+c;
%    w = w.^2+c;
    N = N+1;
    k(~k & abs(z)>4) = N;
%    k(~k & abs(w)>4) = N;
end
k(k==0) = Nmax;
figure
s = pcolor(x, y, mod(k, 2));
colormap([0 0 0;1 1 1])
set(s,'edgecolor','none')
axis([xmin xmax -ymax ymax])
fsize=16;
fsizea=30;
%set(gca,'xtick',[xmin:0.4:xmax],'FontSize',fsize)
%set(gca,'ytick',[-ymax:0.5:ymax],'FontSize',fsize)
xlabel('Real(z)','FontSize',fsizea)
ylabel('Imaginary(z)','FontSize',fsizea)
axis square
title('Julia Set Corresp. to the c you Picked','FontSize',fsizea/2)
rect = [0 rect(2)+10 rect(3) rect(4)];
set(gcf,'Position',rect);
set(gcf,'PaperPositionMode','auto');   

hold on
iterates = 200;		% how many times we iterate the orbit
lag = 15;		% how many iterates are plotted on figure

disp('Pick an initial condition in the Juila Set figure and PRESS ENTER'); pause(1)
disp('Again, a point in the /inside/ of the set will result in more interesting behavior');

[xup,yup] = ginput;
orbit = complex(xup,yup);
disp(sprintf('You chose z = %0.3f + %0.3fi',xup,yup))
disp('Now we will iterate the map f starting with this z')

  for j = 1:iterates
    if mod(j,5) == 0
      disp(sprintf('Iterate %s',int2str(j)))
    end
    if j == 16
      disp('Removing Old Iterates')
    end
    orbit
    orbit = orbit.^2+c;
    h(j)=plot(orbit,'r+','MarkerSize',20);
    if j > lag
      delete(h(j-lag))
    end
    pause(.02)
  end
  disp(''); pause(3)
  disp('What period do you think your point falls in?')
  disp(''); pause(3)
  disp('All points in the interior of the Julia set are attracted to the same sink')
  disp(''); pause(3)
  disp('All points outside are attracted to infinity')
  disp(''); pause(3)
  disp('Press enter to try again')
  pause
  disp('Starting again'); disp(' '); disp(' ');
  for j = iterates-lag+1:iterates
    delete(h(j))
  end
  clear all
  close all
  mandelbrot_points
