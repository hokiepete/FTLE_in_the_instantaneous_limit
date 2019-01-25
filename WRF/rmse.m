close all
clear all
clc
load FTLE_wrf
load correction_term
s1 = s1(:,:,1);
c1 = cor(:,:,1);
sigma(:,:,1)=-s1;
n = length(T);
for i =1:n
    ftle_t = squeeze(sigma(:,:,i));
    sig_true = reshape(ftle_t,[],1);
    sig_approx = reshape(-s1-T(i)*c1,[],1);
    ind = ~isnan(sig_true) & ~isnan(sig_approx) ;
    sig_true = sig_true(ind);
    sig_approx = sig_approx(ind);
    rmse_corrected(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    n=length(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr_corrected(i) = numerator./denominator;

    sig_approx = reshape(-s1,[],1);
    sig_approx = sig_approx(ind);
    rmse_uncorrected(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    corr_uncorrected(i) = numerator./denominator;

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
plot(abs(T),rmse_corrected,'b.-')
plot(abs(T),rmse_uncorrected,'r.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
subplot(122)
hold on
plot(abs(T),corr_corrected,'b.-')
plot(abs(T),corr_uncorrected,'r.-')
legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')


%save wrf_plot_data rmse_corrected rmse_uncorrected time
