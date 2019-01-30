close all
clear all
clc
load FTLE_wrf
load correction3rd
s1 = s1(:,:,1);
l1 = l1(:,:,1);
l2 = l2(:,:,1);
sigma(:,:,1)=-s1;
n = length(T);
for i =1:n
    ftle_t = squeeze(sigma(:,:,i));
    sig_true = reshape(ftle_t,[],1);
    meanf = mean(sig_true);
    sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*l2)),[],1);
    ind = ~isnan(sig_true) & ~isnan(sig_approx) ;
    sig_true = sig_true(ind);
    sig_approx = sig_approx(ind);
    rmse3(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    n=length(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr3(i) = numerator./denominator;
    
    sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1),[],1);
    sig_approx = sig_approx(ind);
    rmse2(i) = sqrt(mean((sig_approx-sig_true).^2));
    mean2(i) = mean(sig_approx);
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr2(i) = numerator./denominator;
    
    sig_approx = reshape(-s1,[],1);
    sig_approx = sig_approx(ind);
    mean1(i) = mean(sig_approx);
    rmse1(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr1(i) = numerator./denominator;

    %cor(i) = corr(sig_approx,sig_true);
end
%time=time(2:27)
%{
start = 1;
stop = 14;
%time = time*60;
time=time(start:stop);
rmse_corrected=rmse_corrected(start:stop);
rmse_uncorrected=rmse_uncorrected(start:stop);
cor_corrected=cor_corrected(start:stop);
cor_uncorrected=cor_uncorrected(start:stop);
%}
figure
subplot(121)
hold on
plot(abs(T),rmse3,'k.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmse1,'r.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
subplot(122)
hold on
plot(abs(T),corr3,'k.-')
plot(abs(T),corr2,'b.-')
plot(abs(T),corr1,'r.-')
legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')


figure
plot(abs(T),meanf,'k.-')
plot(abs(T),mean2,'b.-')
plot(abs(T),mean1,'r.-')
legend('FTLE','-s1-T*corr','-s1','Location','northeast')
ylabel('s^{-1}')
xlabel('|T| s')
title('means')


%
figure
subplot(121)
hold on
plot(abs(T),rmse3,'k.-')
plot(abs(T),rmse2,'b.-')
plot(abs(T),rmse1,'r.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
subplot(122)
hold on
plot(abs(T),corr3,'k.-')
plot(abs(T),corr2,'b.-')
plot(abs(T),corr1,'r.-')
legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%}

%save wrf_plot_data rmse_corrected rmse_uncorrected time
