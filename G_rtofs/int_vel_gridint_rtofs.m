close all
clear all
clc
load rtofs_vel
xwant = 400
ywant= 320
time= linspace(0,24*3600,25)

[xx,yy,tt]=meshgrid(x,y,time);

t0 = 24*3600;
tf = 0*3600;
tdim = 97;
t_want = linspace(t0,tf,tdim);

P = [2,1,3];
u = permute(u,P);
v = permute(v,P);
xx = permute(xx,P);
yy = permute(yy,P);
tt = permute(tt,P);
U = griddedInterpolant(xx,yy,tt,u,'cubic','none')
V = griddedInterpolant(xx,yy,tt,v,'cubic','none')

[xx,yy]=meshgrid(x,y);
fx = NaN(ywant,xwant,tdim);
fy = NaN(ywant,xwant,tdim);
warning('off','all')
TSPAN = [tf t0]; % Solve from t=1 to t=5
opts=odeset('event',@eventfun_gridint_rtofs,'RelTol',1e-14)%,'AbsTol',1e-14);

tic
Y0=[xx(160,100),yy(160,100)];
%[T Y] = ode45(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
[T Y] = ode113(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
%{
for i =1:ywant
    for j=1:xwant
        sprintf('%03d, %03d',i,j)
        tic
        Y0=[xx(i,j),yy(i,j)];
        %[T Y] = ode45(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
        [T Y] = ode113(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
        if sum(~isnan(Y(:,1)))<2
            fx(i,j,1) = Y(1,1);
        else
            fx(i,j,:) = interp1(T,squeeze(Y(:,1)),t_want,'spline',nan);
        end
        if sum(~isnan(Y(:,2)))<2
            fy(i,j,1) = Y(1,2);
        else
            fy(i,j,:) = interp1(T,squeeze(Y(:,2)),t_want,'spline',nan);
        end

        toc
    end
end
time = t_want;
save flow_map_gridint_rtofs fx fy xx yy time
%}