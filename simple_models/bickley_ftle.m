clear all
close all
clc
tlen = 101
t0 = 0
tf = -0.2
ydim = 101;
xdim = ceil(2.5*ydim)
y=linspace(-4000,4000,ydim); %km
dy = y(2)-y(1);
x=linspace(0,20000,xdim); %km
dx = x(2)-x(1);
[x,y]=meshgrid(x,y);
twant = linspace(t0,tf,tlen);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
for i = 1:ydim
    i
    for j = 1:xdim
            y0=[x(i,j),y(i,j)];
            [t,yout] = ode45(@bickley_int,[t0,tf],y0,options);
            fx(i,j,:) = interp1(t,yout(:,1),twant,'spline');
            fy(i,j,:) = interp1(t,yout(:,2),twant,'spline');
    end
end
sigma = zeros([ydim,xdim,tlen]);
size(sigma)
for t = 1:tlen
    
    %need to skip t0 to avoid division by zero
    %set FTLE to zero instead
    if twant(t) == t0
        continue
    else
        [dfxdx,dfxdy] = gradient(squeeze(fx(:,:,t)),dx,dy);
        [dfydx,dfydy] = gradient(squeeze(fy(:,:,t)),dx,dy);
        for i = 1:ydim
            for j = 1:xdim
                    gradF = [dfxdx(i,j),dfxdy(i,j);
                             dfydx(i,j),dfydy(i,j)];
                    C = gradF'*gradF;
                    lambda = max(eig(C));
                    sigma(i,j,t) = 1/(2*abs(twant(t)-t0))*log(lambda);
            end
        end
    end
end
size(sigma)
time=twant;
save('bickley_sigma_big.mat', 'sigma','time');
