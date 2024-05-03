% Chapter 5 - Fractals and Multifractals.
% Copyright Birkhauser 2004. Stephen Lynch.
% Modified by CD 2007

% Plot the Koch curve up to stage k (Figure 5.2)

clear all
close all

for k = 1:6
mmax=4^k;
x=zeros(1,mmax);
y=zeros(1,mmax);
h=3^(-k);
x(1)=0;y(1)=0;
angle(1)=0;angle(2)=pi/3;angle(3)=-pi/3;angle(4)=0;
for a=1:mmax
    m=a-1;ang=0;
    for b=1:k
        segment(b)=mod(m,4);
        m=floor(m/4);
        r=segment(b)+1;
        ang=ang+angle(r);
    end
    x(a+1)=x(a)+h*cos(ang);
    y(a+1)=y(a)+h*sin(ang);
end

plot(x,y,'b');
axis equal
pause
end

figure
dim = boxcount([x;y]',9,1);
