
function int_vel()
close all
clear all
clc
load wrf_vel_data
xwant = 405
ywant= 325
t_want = linspace(24,22,11);
yy=linspace(0,972000,ywant);
xx=linspace(0,1212000,xwant);
[xx,yy]=meshgrid(xx,yy);
TSPAN = [24 22]; % Solve from t=1 to t=5
opts=odeset('event',@eventfun,'RelTol',1e-8,'AbsTol',1e-8);

for i =1:ywant
    for j=1:xwant
        sprintf('%03d, %03d',i,j)
        Y0=[xx(i,j),yy(i,j)];
        [T Y] = ode45(@odefun, TSPAN, Y0, opts, x, y, time,u,v); % Solve ODE
        fx(i,j,:) = interp1(T,squeeze(Y(:,1)),t_want,'spline',nan);
        fy(i,j,:) = interp1(T,squeeze(Y(:,2)),t_want,'spline',nan);
    end
end
time = t_want;
save flow_map fx fy xx yy time
end


function dydt = odefun(t,Y,x,y,time,u,v)
    dydt = zeros(2,1);
    dydt(1) = interp3(x,y,time,u,Y(1),Y(2),t,'spline');
    dydt(2) = interp3(x,y,time,v,Y(1),Y(2),t,'spline');
end

function [value,isterminal,direction]=eventfun(t,Y,x,y,time,u,v)
    value=[Y(1),Y(1)-1212000,Y(2),Y(2)-972000];
    isterminal=[1,1,1,1];
    direction=[0,0,0,0];
end
