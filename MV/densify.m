close all
clear all
clc
load mv_vel_data
uu = u;
vv = v;
clear u v

t_want = [0:200:24*3600];
for i = 1:length(y)
    i
    for j = 1:length(x);
        u(i,j,:) = interp1(time,squeeze(uu(i,j,:)),t_want,'spline');
        v(i,j,:) = interp1(time,squeeze(vv(i,j,:)),t_want,'spline');
    end
end
time = t_want;
save mv_vel_data_dense x y time u v