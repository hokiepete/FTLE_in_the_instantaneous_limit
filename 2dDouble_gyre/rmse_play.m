close all
clear all
clc
load dg_sigma_big_short
load error_comparison
s1 = s1(:,:,2);
l1 = l1(:,:,2);
a1 = a1(:,:,2);
a2 = a2(:,:,2);

a2(abs(a2)>7) = nan;
sigma(:,:,1)=-s1;
T=time;
n = length(T);
for i =1:n
    ftle_t = squeeze(sigma(:,:,i));
    sig_true = reshape(ftle_t,[],1);
    sig_approx2a = reshape(-s1-T(i).*(-s1.^2+0.5.*l1-T(i)),[],1);
    sig_approx2b = reshape(-s1-T(i).*(-s1.^2+0.5.*l1-T(i).*(4/3.*s1.^3)),[],1);
    sig_approx1 = reshape(-s1-T(i).*(-s1.^2+0.5.*l1-T(i).*(-s1.*l1)),[],1);
    s1_approx = reshape(-s1-T(i).*(-s1.^2+0.5.*l1-T(i).*(0.25.*a2)),[],1);
    sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a2)),[],1);
    sig_cor = reshape(-s1-T(i)*(-s1.^2+0.5*l1),[],1);   
    
    %sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a1)),[],1);
    %sig_approx2 = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a2)),[],1);   
    
    ind = ~isnan(sig_cor) & ~isnan(sig_approx) & ~isnan(sig_true) & ~isnan(sig_approx1) & ~isnan(sig_approx2a) & ~isnan(sig_approx2b) & ~isinf(sig_cor) & ~isinf(sig_approx) & ~isinf(sig_true) & ~isinf(sig_approx1) & ~isinf(sig_approx2a) & ~isinf(sig_approx2b);
    sig_true = sig_true(ind);
    
    sig_approx2a = sig_approx2a(ind);
    rmse2a(i) = sqrt(mean((sig_approx2a-sig_true).^2));
    
    sig_approx2b = sig_approx2b(ind);
    rmse2b(i) = sqrt(mean((sig_approx2b-sig_true).^2));
    
    sig_approx1 = sig_approx1(ind);
    rmse1(i) = sqrt(mean((sig_approx1-sig_true).^2));
    
    sig_approx = sig_approx(ind);
    rmse(i) = sqrt(mean((sig_approx-sig_true).^2));

    s1_approx = s1_approx(ind);
    rmse0(i) = sqrt(mean((s1_approx-sig_true).^2));

    sig_cor = sig_cor(ind);
    rmsecor(i) = sqrt(mean((sig_cor-sig_true).^2));

end

figure
hold on
plot(abs(T),rmse0,'r.-')
plot(abs(T),rmse1,'b.-')
plot(abs(T),rmse2a,'k.-')
plot(abs(T),rmse2b,'g.-')
plot(abs(T),rmse,'c.-')
plot(abs(T),rmsecor,'m.-')
%legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
%legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0-d.^2/m')

%
figure
hold on
plot(abs(T),rmse0,'r.-')
plot(abs(T),rmse1,'b.-')
plot(abs(T),rmse2a,'k.-')
plot(abs(T),rmse2b,'g.-')
plot(abs(T),rmse,'c.-')
plot(abs(T),rmsecor,'m.-')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
%}

%{
for i =1:n
    ftle_t = squeeze(sigma(:,:,i));
    sig_true = reshape(ftle_t,[],1);
    sig_approx2a = reshape(-s1-T(i).*(-s1.^2+0.5.*l1-T(i)),[],1);
    sig_approx2b = reshape(-s1-T(i).*(-s1.^2+0.5.*l1-T(i).*(4/3.*s1.^3-s1.*l1+0.25.*a2)),[],1);
    sig_approx1 = reshape(-s1-T(i).*(-s1.^2+0.5.*l1),[],1);
    s1_approx = reshape(-s1,[],1);
    %sig_approx = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a1)),[],1);
    %sig_approx2 = reshape(-s1-T(i)*(-s1.^2+0.5*l1+T(i)*(4/3*s1.^3-s1.*l1+0.25*a2)),[],1);   
    
    ind = ~isnan(sig_true) & ~isnan(sig_approx1) & ~isnan(sig_approx2a) & ~isnan(sig_approx2b) & ~isinf(sig_true) & ~isinf(sig_approx1) & ~isinf(sig_approx2a) & ~isinf(sig_approx2b);
    sig_true = sig_true(ind);
    
    sig_approx2a = sig_approx2a(ind);
    rmse2a(i) = sqrt(mean((sig_approx2a-sig_true).^2));
    
    sig_approx2b = sig_approx2b(ind);
    rmse2b(i) = sqrt(mean((sig_approx2b-sig_true).^2));
    
    sig_approx1 = sig_approx1(ind);
    rmse1(i) = sqrt(mean((sig_approx1-sig_true).^2));

    s1_approx = s1_approx(ind);
    rmse0(i) = sqrt(mean((s1_approx-sig_true).^2));
   
end

figure
subplot(121)
hold on
plot(abs(T),rmse0,'r.-')
plot(abs(T),rmse1,'b.-')
plot(abs(T),rmse2a,'k.-')
%legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')

subplot(122)
hold on
plot(abs(T),rmse0,'r.-')
plot(abs(T),rmse1,'b.-')
plot(abs(T),rmse2b,'k.-')
%legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0-d.^2/m')

%
figure
subplot(121)
hold on
plot(abs(T),rmse0,'r.-')
plot(abs(T),rmse1,'b.-')
plot(abs(T),rmse2a,'k.-')
%legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
subplot(122)
hold on
plot(abs(T),rmse0,'r.-')
plot(abs(T),rmse1,'b.-')
plot(abs(T),rmse2b,'k.-')
%legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
%}

%save wrf_plot_data rmse_corrected rmse_uncorrected time
%}