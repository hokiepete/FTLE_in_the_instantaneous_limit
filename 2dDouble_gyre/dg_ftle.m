clear all
close all
clc

%Set up the domain
lenx = 301
leny = 151
tlen = 121
t0 = 0
tf = -0.8
x = linspace(0,2,lenx);
dx=x(2)-x(1);
y = linspace(0,1,leny);
dy=y(2)-y(1);
time = linspace(t0,tf,tlen);
[x,y] = meshgrid(x,y);

%Calculate the flow map
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
for i = 1:leny
    for j = 1:lenx
            y0=[x(i,j),y(i,j)];
            [t,yout] = ode113(@dg_int,[t0,tf],y0,options);
            fx(i,j,:) = interp1(t,yout(:,1),time,'spline');
            fy(i,j,:) = interp1(t,yout(:,2),time,'spline');
    end
end

%Compute the FTLE
sigma = zeros([leny,lenx,tlen]);
for t = 1:tlen
    %need to skip t0 to avoid division by zero
    %set FTLE to zero instead
    if time(t) == t0
        continue
    else
        [dfxdx,dfxdy] = gradient(squeeze(fx(:,:,t)),dx,dy);
        [dfydx,dfydy] = gradient(squeeze(fy(:,:,t)),dx,dy);
        for i = 1:leny
            for j = 1:lenx
                    gradF = [dfxdx(i,j),dfxdy(i,j);
                             dfydx(i,j),dfydy(i,j)];
                    C = gradF'*gradF;
                    lambda = max(eig(C));
                    sigma(i,j,t) = 1/(2*abs(time(t)-t0))*log(lambda);
            end
        end
    end
end

save('dg_sigma.mat', 'sigma','time');
