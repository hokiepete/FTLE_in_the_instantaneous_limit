%function int_vel_gridint_v2()
clear all
close all
clc
load mv_vel_data_dense

t0 = 24*3600;
tf = 0*3600;
tdim = 241;
ydim = length(y);
xdim = length(x);
t_want = linspace(tf,t0,tdim);
[xx,yy,tt] = meshgrid(x,y,time);
[x,y] = meshgrid(x,y);

%
P = [2,1,3];
u = permute(u,P);
v = permute(v,P);
xx = permute(xx,P);
yy = permute(yy,P);
tt = permute(tt,P);
%}
u(u==999)=nan;
v(v==999)=nan;
U = griddedInterpolant(xx,yy,tt,u,'cubic','none')
V = griddedInterpolant(xx,yy,tt,v,'cubic','none')

%
fx = NaN(ydim,xdim,tdim);
fy = NaN(ydim,xdim,tdim);
TSPAN = [t0,tf]; % Solve from t=1 to t=5
warning('off','all')
opts=odeset('event',@eventfun_gridint_mv,'RelTol',1e-14)%,'AbsTol',1e-14);
for i =1:ydim
    for j=1:xdim
        sprintf('%03d, %03d',i,j)
        tic
        Y0=[x(i,j),y(i,j)];
        %[T Y] = ode45(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
        [T Y] = ode113(@odefun_gridint, TSPAN, Y0, opts, U,V); % Solve ODE
        if sum(~isnan(Y(:,1)))<2
            fx(i,j,end) = Y(1,1);
        else
            fx(i,j,:) = interp1(T,squeeze(Y(:,1)),t_want,'spline',nan);
        end
        if sum(~isnan(Y(:,2)))<2
            fy(i,j,end) = Y(1,2);
        else
            fy(i,j,:) = interp1(T,squeeze(Y(:,2)),t_want,'spline',nan);
        end
        toc
    end
end
time = t_want;
save flow_map_gridint_v2_mv_sc fx fy x y time

%{
end


function [value,isterminal,direction]=eventfun_gridint_mv(t,Y,U,V)
    value=[isnan(Y(1))-1,isnan(Y(2))-1];
    isterminal=[1,1];
    direction=[0,0];
end

%
function [value,isterminal,direction]=eventfun_gridint_mv(t,Y,U,V)
    value=[Y(1),Y(1)-52400,Y(2),Y(2)-47800,U(Y(1),Y(2),t)-2.270778503417969,V(Y(1),Y(2),t)-2.500529632568359];
    isterminal=[1,1,1,1,1,1];
    direction=[0,0,0,0,+1,+1];
end
%}
%}