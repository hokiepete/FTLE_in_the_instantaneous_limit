y=0:0.01:0.5;
y=y(:);

k=0;
for T=0.001:0.001:0.2,
k=k+1;
T=-T;
sig_true=-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3));
s1=-(1+3*y.^2);
lam1=2 + 18*y.^2 + 24*y.^4 ;
lam2=-(8/3 + 56*y.^2 + 192*y.^4 + 160*y.^6);

sig_0 = -s1;
sig_1 = -(-s1.^2 + 1/2*lam1) * T ;
sig_2 = -(4/3*s1.^3 - s1.*lam1 + 1/4*lam2) * T^2 ;
err0(k)=sqrt(immse(sig_true,sig_0));
rho0(k)=corr(      sig_true,sig_0);
err1(k)=sqrt(immse(sig_true,sig_0+sig_1));
err2(k)=sqrt(immse(sig_true,sig_0+sig_1+sig_2));
rho1(k)=corr(      sig_true,sig_0+sig_1);
t(k)=T;
end

plot(abs(t),err0)
hold on
plot(abs(t),err1,'m')
plot(abs(t),err2,'k')
xlabel('|T|');
ylabel('Root mean-squared error');

figure
plot(abs(t),rho0)
hold on
plot(abs(t),rho1,'m')
xlabel('|T|');
ylabel('Pearson correlation coefficient');