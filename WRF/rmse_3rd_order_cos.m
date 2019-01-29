close all
clear all
clc

n = 101;
x = linspace(0,24,n);
true = 1./(x+1);
approx1 = 1;
aerror1 = abs(approx1-true);
rerror1 = abs(approx1-true)./abs(true);

approx2 = 1-x;
aerror2 = abs(approx2-true);
rerror2 = abs(approx2-true)./abs(true);

approx3 = 1-x+x.^2;
aerror3 = abs(approx3-true);
rerror3 = abs(approx3-true)./abs(true);

figure
hold on
plot(x,true,'g.-')
plot(x,approx3,'k.-')
plot(x,approx2,'b.-')
plot(x,approx1,'r.-')

legend('true','3rd order','2nd order','1st order','Location','northwest')


figure
subplot(121)
hold on
plot(x,aerror3,'k.-')
plot(x,aerror2,'b.-')
plot(x,aerror1,'r.-')
legend('3rd order','2nd order','1st order','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
subplot(122)
hold on
plot(x,rerror3,'k.-')
plot(x,rerror2,'b.-')
plot(x,rerror1,'r.-')
legend('3rd order','2nd order','1st order','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')

%

%{
figure
subplot(121)
hold on
plot(x,aerror3,'k.-')
plot(x,aerror2,'b.-')
plot(x,aerror1,'r.-')
legend('3rd order','2nd order','1st order','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
subplot(122)
hold on
plot(x,rerror3,'k.-')
plot(x,rerror2,'b.-')
plot(x,rerror1,'r.-')
legend('3rd order','2nd order','1st order','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%}

%save wrf_plot_data aerror_rerrorected aerror_unrerrorected time
%{
true = 1./(x+1).*cos(x);
approx = 1;
aerror1 = abs(approx-true);
rerror1 = abs(approx-true)./abs(true);

approx = 1-x;
aerror2 = abs(approx-true);
rerror2 = abs(approx-true)./abs(true);

approx = 1-x+0.5*x.^2;
aerror3 = abs(approx-true);
rerror3 = abs(approx-true)./abs(true);

%}