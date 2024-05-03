function [LE] = LEcalc(f,Df,param,v0,type)

if nargin<5, type = 'map'; end

m = length(v0);
N = 10^3;
LE = zeros(m,1);
v = v0;

if strcmp(type,'map')
    
    Z = Df(v,param);
    [q,r] = qr(Z);
    LE = LE + log(abs(diag(r)))/N;
    for i=2:N
        v = f(v,param);
        Z = Df(v,param)*q;
        [q,r] = qr(Z);
        LE = LE + log(abs(diag(r)))/N;
    end
    %LE = LE/N;
    
elseif strcmp(type,'ode')
    
    tstep = 0.01;
    q = eye(m);
    Jt = @(x,A) A*x; % ode for DF1, 'A' will be current Df
    
    for i=1:N
        DF1 = eye(m); % initialize derivative of time-1 map
        for j = 1:1/tstep
            % one timestep using RK4
            DF1 = stepit(Jt,DF1,Df(v,param),tstep);	% 3 x 3 matrix
            v = stepit(f,v,param,tstep);
        end
        disp(DF1)
        Z = DF1*q;
        [q,r] = qr(Z);
        LE = LE + log(abs(diag(r)))/N;
    end
    %LE = LE/N;
else
    
    disp('Unrecognized type. LE = 0-vector returned.')
   
end

end
