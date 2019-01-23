function int_vel_gridint_v2()
close all
clear all
clc
load wrf_vel_data
xwant = 51%405
ywant= 51%325
t_want = linspace(24*3600,22*3600,121);
yy=linspace(0,972000,ywant);
xx=linspace(0,1212000,xwant);
P = [2,1,3];
u = permute(u,P);
v = permute(v,P);
x = permute(x,P);
y = permute(y,P);
time = permute(time,P);
U = griddedInterpolant(x,y,time,u,'spline','spline')
V = griddedInterpolant(x,y,time,v,'spline','spline')

[xx,yy]=meshgrid(xx,yy);
fx = NaN(ywant,xwant,121);
fy = NaN(ywant,xwant,121);
TSPAN = [24*3600 22*3600]; % Solve from t=1 to t=5
opts=odeset('event',@eventfun_gridint,'RelTol',1e-8)%,'AbsTol',1e-8);
for i =1:ywant
    for j=1:xwant
        sprintf('%03d, %03d',i,j)
        tic
        Y0=[xx(i,j),yy(i,j)];
        [T Y] = ode45(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
        fx(i,j,:) = interp1(T,squeeze(Y(:,1)),t_want,'spline',nan);
        fy(i,j,:) = interp1(T,squeeze(Y(:,2)),t_want,'spline',nan);
        toc
    end
end
time = t_want;
save flow_map_gridint_v2 fx fy xx yy time
end

function dydt = odefun_gridint(t,Y,U,V)
    %[Y(1),Y(2),t]
    %[Y(1),Y(1)-1212000,Y(2),Y(2)-972000]
    dydt = zeros(2,1);
    dydt(1) = U(Y(1),Y(2),t);%interp3(x,y,time,u,Y(1),Y(2),t,'spline');
    dydt(2) = V(Y(1),Y(2),t);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');
end

function [value,isterminal,direction]=eventfun_gridint(t,Y,U,V)
    value=[Y(1)+0.0000001,Y(1)-1212000.0000001,Y(2)+0.0000001,Y(2)-972000.0000001];
    isterminal=[1,1,1,1];
    direction=[0,0,0,0];
end
