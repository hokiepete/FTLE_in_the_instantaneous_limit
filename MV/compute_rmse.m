close all
clear all
clc

load error_comparison_terms

load FTLE_mv_short
%time_ftle = T
[ydim,xdim,tdim]=size(sigma)
j = 13
s1 = s1(:,:,j);
l1 = l1(:,:,j);
a1 = a1(:,:,j);
a2 = a2(:,:,j);

i=0;
s1 = s1((i+1):end-i,(i+1):end-i,1);
l1 = l1((i+1):end-i,(i+1):end-i,1);
a1 = a1((i+1):end-i,(i+1):end-i,1);
a2 = a2((i+1):end-i,(i+1):end-i,1);
sigma = sigma((i+1):end-i,(i+1):end-i,:);

sigma(:,:,1)=-s1;
n = length(T);%time_ftle);
for i =1:n
    %T = time_ftle(end-i+1);
    ftle_t = squeeze(sigma(:,:,i));
    sig_true = reshape(ftle_t,[],1);
    sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a1)),[],1);
    sig_approx2 = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a2)),[],1);
    
    ind = ~isnan(sig_true) & ~isnan(sig_approx) & ~isnan(sig_approx2);
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

figure
subplot(121)
hold on
plot(-T,rmse1,'r.-')
plot(-T,rmse2,'b.-')
plot(-T,rmsea1,'k.-')
legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')

subplot(122)
hold on
plot(-T,rmse1,'r.-')
plot(-T,rmse2,'b.-')
plot(-T,rmsea2,'k.-')
%xlim([0,250])
legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0-d.^2/m')

%
figure
subplot(121)
hold on
plot(-T,rmse1,'r.-')
plot(-T,rmse2,'b.-')
plot(T,rmsea1,'k.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
subplot(122)
hold on
%plot(-T,rmse1,'r.-')
plot(-T,rmse2,'b.-')
plot(T,rmsea2,'k.-')
legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
%}
rmse3 = rmsea2;
time = T;
save wrf_plot_data rmse1 rmse2 rmse3 time


