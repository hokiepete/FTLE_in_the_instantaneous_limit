close all
clear all
clc
load abc_sigma_big%_big_ex_short
%load abc_sigma_shorttime
load abc_error_comparison
%

s1 = s1(:,:,:,2);
l1 = l1(:,:,:,2);
a1 = a1(:,:,:,2);
%a1(abs(a1)>10)=nan;
%a1(abs(a1)>1e1)=nan;
figure
sigma(1,:,:,:)=-s1;
T=twant;
n = length(T);
for i =1:n
    ftle_t = squeeze(sigma(i,:,:,:));
    sig_true = reshape(ftle_t,[],1);
    %sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a1)),[],1);
    %sig_approx2 = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a2)),[],1);
    pt(i) = T(i);
    pt2(i) = T(i).^2;
    
    sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1)-T(i).^2*(4/3*s1.^3-s1.*l1+0.25*a1),[],1);
    %sig_approx2 = reshape(-s1-T(i)*(-s1.^2+0.5*l1)-T(i).^2*(4/3*s1.^3-s1.*l1+0.25*a2),[],1);
    
    
    ind = ~isnan(sig_true) & ~isnan(sig_approx) & ~isinf(sig_true) & ~isinf(sig_approx);
    
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
    
    sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1),[],1);
    sig_approx = sig_approx(ind);
    rmse2(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr2(i) = numerator./denominator;
    
    sig_approx = reshape(-s1,[],1);
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
hold on
%plot(abs(T),rmse1,'r.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmsea1,'k.-')
%plot(abs(T),abs(T).^1,'r--')
%plot(abs(T),m2*abs(T).^2,'b--')
%plot(abs(T),m3*abs(T).^3,'k--')

%xlim([0.,0.5])
legend('-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')

%{
figure
subplot(121)
hold on
plot(abs(T),rmse1,'r.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmsea1,'k.-')
legend('-s1-T*corr','-s1','Location','southeast')
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
plot(abs(T),rmsea1,'k.')
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
rmse3=rmsea1;
time = T;
%save abc_plot_data_shortime rmse1 rmse2 rmse3 time
save abc_plot_data_big rmse1 rmse2 rmse3 time
%}
%}
%}