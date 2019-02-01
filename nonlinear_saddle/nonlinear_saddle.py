# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 22:20:26 2019

@author: pnola
"""

from numpy import n

y=0:0.01:0.5;
y=y(:);

k=0;
for T=0.01:0.01:0.8,
k=k+1;
T=-T;
sig_0 = (1+3*y.^2);
sig_true=-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3));
err(k)=immse(sig_true,sig_0);
rho(k)=corr(sig_true,sig_0);
t(k)=T;
end

k=0;
for T=0.01:0.01:0.8,
k=k+1;
T=-T;
sig_0 = (1+3*y.^2);
sig_1 = -3*T*y.^2.*(1+3*y.^2);
sig_true=-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3));
err2(k)=immse(sig_true,sig_0+sig_1);
rho2(k)=corr(sig_true,sig_0+sig_1);
t(k)=T;
end

plot(-t,rho2,'m')
hold on
plot(-t,rho)
figure
plot(-t,err)
hold on
plot(-t,err2,'m')
figure(2);
xlabel('|T|');
ylabel('Mean-squared error');

figure(1)
xlabel('|T|');
ylabel('Pearson correlation coefficient');