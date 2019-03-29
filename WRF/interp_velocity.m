%{
function interp_velocity()
close all
clear all
clc
load wrf_vel_data
xwant = 405
ywant= 325
t_want = 24*3600%linspace(24*3600,22*3600,121);
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

[xx,yy,t_want]=meshgrid(xx,yy,t_want);
%{
xx = permute(xx,P);
yy = permute(yy,P);
t_want = permute(t_want,P);
%}

uu = U(xx,yy,t_want);%interp3(x,y,time,u,Y(1),Y(2),t,'spline');
vv = V(xx,yy,t_want);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');
time = t_want;
%{
uu = permute(uu,P);
vv = permute(vv,P);
xx = permute(xx,P);
yy = permute(yy,P);
time = permute(time,P);
%}
save high_res_vel uu vv xx yy time
end
%}

function interp_velocity()
close all
clear all
clc
load wrf_vel_data
xwant = 405
ywant= 325
t0 = 13*3600;
tf = 11*3600;

t_want = linspace(t0,tf,3);
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

[xx,yy,t_want]=meshgrid(xx,yy,t_want);
xx = permute(xx,P);
yy = permute(yy,P);
t_want = permute(t_want,P);


uu = U(xx,yy,t_want);%interp3(x,y,time,u,Y(1),Y(2),t,'spline');
vv = V(xx,yy,t_want);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');
time = t_want;

uu = permute(uu,P);
vv = permute(vv,P);
xx = permute(xx,P);
yy = permute(yy,P);
%time = permute(time,P);
xx=xx(:,:,1);
yy=yy(:,:,1);
time=squeeze(time(1,1,:));
save high_res_vel uu vv xx yy time
end
%}