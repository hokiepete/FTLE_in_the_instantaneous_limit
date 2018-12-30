clear all
close all
clc
len = 96
t0 = 0
tf = 2*pi
x = linspace(0,2*pi,len);
dx=x(2)-x(1);
twant = fliplr(linspace(t0,tf,len));
load output
t=2
figure
vol3d_v2('cdata',squeeze(sigma(t,:,:,:)))
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
%input('ready?')

load s1
figure
vol3d_v2('cdata',-s1)
xlabel('x')
ylabel('y')
zlabel('z')
colorbar


%}
