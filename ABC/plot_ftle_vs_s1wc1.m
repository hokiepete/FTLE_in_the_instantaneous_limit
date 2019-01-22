close all
clear all
clc
tlen = 96
t0 = 0
tf = -2.5
time = linspace(t0,tf,tlen);
load sigma_big
load eulerian_rmse_data
c1 = correction;
ftle = sigma;
ftle(1,:,:,:)=-s1;
n = length(time);
for i =1:4
    T=time(i);
    ftle_t = squeeze(ftle(i,:,:,:));
    sig_true = ftle_t;
    sig_approx = -s1-T*c1;
    figure
    subplot(211)
    vol3d_v2('cdata',sig_true)
    colorbar
    subplot(212)
    vol3d_v2('cdata',sig_approx)
    colorbar
end