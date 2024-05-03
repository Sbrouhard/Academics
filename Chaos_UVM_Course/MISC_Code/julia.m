% Copyright Birkhauser 2004. Stephen Lynch.
% Commented/modified by Chris Danforth 3/16/10
% Plot the Julia set J(0,0.9)

clear all
close all

k=20; niter=2^k;			% 2^k total iterations
x=zeros(1,niter); y=zeros(1,niter);	% filling empty arrays
x1=zeros(1,niter); y1=zeros(1,niter);	
%a=0; b=.9;				% real and imaginary parts of c
a=0; b=1.1;				% real and imaginary parts of c
x(1)=real(0.5+sqrt(0.25-(a+i*b)));	% IC
y(1)=imag(0.5+sqrt(0.25-(a+i*b)));	% IC
% isunstable=2*abs(x(1)+i*y(1));	% Check that the point is unstable, don't need this now.

for n=1:niter				% iterate 
    x1=x(n); y1=y(n);
    u=sqrt((x1-a)^2+(y1-b)^2)/2; v=(x1-a)/2;
    u1=sqrt(u+v); v1=sqrt(u-v);
    x(n+1)=u1; y(n+1)=v1;
    if y(n)<b 
      y(n+1)=-y(n+1);
    end
    if rand < .5			% flip a coin
      x(n+1)=-u1; y(n+1)=-y(n+1);
    end
end

plot(x,y,'k.','MarkerSize',2)

%% trick out the figure before printing
fsize=15; fsizea=30;			% font size of axis and title labels
xlim([-1.6 1.6]); ylim([-1.6 1.6]);
axis square
xlabel('Real(z)','FontSize',fsizea)
ylabel('Imaginary(z)','FontSize',fsizea)
set(gca,'xtick',[-1.6:0.4:1.6],'FontSize',fsize)
set(gca,'ytick',[-1.2:0.4:1.2],'FontSize',fsize)
title(sprintf('Julia Set for c = %s + %s i',num2str(a,'%1.1d'),num2str(b,'%1.1d')),'FontSize',fsizea)
rect = [0,0,260*2.0,260*1.7];
%rect = [0,0,560*2.0,560*1.7];
set(gcf,'Position',rect);
set(gcf,'PaperPositionMode','auto');   
%print('-djpeg','picture_julia')    
