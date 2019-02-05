close all
clear all
clc
%load dg_sigma_big_short
load dg_sigma_big
load analytic_lambda_terms
%
%
figure
subplot(121)
surface(-lambda_0,'edgecolor','none')
colorbar
subplot(122)
surface(sigma(:,:,3),'edgecolor','none')
colorbar
sigma(:,:,1)=-lambda_0;
T=time;
n = length(T);

for i =1:n
    ftle_t = squeeze(sigma(:,:,i));
    sig_true = reshape(ftle_t,[],1);
    sig_approx = reshape(-lambda_0-T(i)*(-lambda_0.^2+0.5*lambda_1+T(i)*(4/3*lambda_0.^3-lambda_0.*lambda_1+0.25*lambda_2_first)),[],1);
    sig_approx2 = reshape(-lambda_0-T(i)*(-lambda_0.^2+0.5*lambda_1+T(i)*(4/3*lambda_0.^3-lambda_0.*lambda_1+0.25*lambda_2_second)),[],1);
    
    %sig_approx = reshape(-lambda_0-T(i)*(-lambda_0.^2+0.5*lambda_1)-T(i).^2*(4/3*lambda_0.^3-lambda_0.*lambda_1+0.25*lambda_2_first),[],1);
    %sig_approx2 = reshape(-lambda_0-T(i)*(-lambda_0.^2+0.5*lambda_1)-T(i).^2*(4/3*lambda_0.^3-lambda_0.*lambda_1+0.25*lambda_2_second),[],1);
    
    
    ind = ~isnan(sig_true) & ~isnan(sig_approx) & ~isnan(sig_approx2) & ~isinf(sig_true) & ~isinf(sig_approx) & ~isinf(sig_approx2);
    sig_true = sig_true(ind);
    sig_approx = sig_approx(ind);
    rmsea1(i) = sqrt(mean((sig_approx-sig_true).^2));
    
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    n=length(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corra1(i) = numerator./denominator;
    
    sig_approx = sig_approx2(ind);
    rmsea2(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    n=length(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corra2(i) = numerator./denominator;
    
    sig_approx = reshape(-lambda_0-T(i)*(-lambda_0.^2+0.5*lambda_1),[],1);
    sig_approx = sig_approx(ind);
    rmse2(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr2(i) = numerator./denominator;
    
    sig_approx = reshape(-lambda_0,[],1);
    sig_approx = sig_approx(ind);
    rmse1(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr1(i) = numerator./denominator;

   
end
%{
T=T(2:end);
rmse1=rmse1(2:end);
rmse2=rmse2(2:end);
rmsea2=rmsea2(2:end);
rmsea1=rmsea1(2:end);
%}
m2=1
m3=1
figure
subplot(121)
hold on
plot(abs(T),rmse1,'r.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmsea1,'k.-')
%{
plot(abs(T),abs(T).^1,'r--')
plot(abs(T),abs(T).^2,'b--')
plot(abs(T),abs(T).^3,'k--')
xlim([0.,0.5])
%}
legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')
subplot(122)
hold on
plot(abs(T),rmse1,'r.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmsea2,'k.-')
%{
plot(abs(T),abs(T).^1,'r--')
plot(abs(T),abs(T).^2,'b--')
plot(abs(T),abs(T).^3,'k--')
xlim([0,0.5])
%}
legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0-d.^2/m')

%{
figure
subplot(121)
hold on
plot(abs(T),rmse1,'r.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmsea1,'k.-')
legend('-lambda_0-T*corr','-lambda_0','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
subplot(122)
%
T=T(2:end)
rmse1=rmse1(2:end)
rmse2=rmse2(2:end)
rmsea2=rmsea2(2:end)
%}

figure
hold on
plot(abs(T),rmse1,'r.')
plot(abs(T),rmse2,'b.')
plot(abs(T),rmsea2,'k.')
plot(abs(T),abs(T).^1,'r--')
plot(abs(T),abs(T).^2,'b--')
plot(abs(T),abs(T).^3,'k--')

legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
axis equal
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%}

%save wrf_plot_data rmse_corrected rmse_uncorrected time
%}
%}
%}