clear all
close all
clc
y=-0.5:0.005:0.5;
dy = y(2)-y(1);
dx = dy;%y(2)-y(1);
dt=1;
ydim = length(y);
xdim=ydim;
[x,y]=meshgrid(y,y);
%y=y(:);
R = [0,-1;1,0];
for i =1:ydim
    for j = 1:xdim
        i,j
        S = [1,0;0,-(1+3*y(i,j)^2)];
        B = [2,0;0,(2+18*y(i,j)^2+24*y(i,j)^4)];
        Q = [8/3,0;0,-(8/3+56*y(i,j)^2+192*y(i,j)^4+160*y(i,j)^6)];
        [V,D] = eig(S);
        if ~issorted(diag(D))
            [D,I] = sort(diag(D));
            V = V(:, I);
        end
        s1(i,j) = D(1,1);
        X0 = V(:,1);
        l1(i,j) = X0'*B*X0;
                
        X1 = -((B-l1(i,j)*eye(size(B)))*X0)\(S-s1(i,j)*eye(size(B)));
        if sum(X1)~=0
            X1=X1'/norm(X1);
        else
            X1=X1';
        end
        l2_a(i,j) = X0'*Q*X0 + X0'*B*X1 - X0'*S*X1;
        m = X0'*R'*(S-s1(i,j)*eye(size(S)))*R*X0;
        d = X0'*R'*B*X0;
        l2_b(i,j) = X0'*Q*X0-d.^2/m;
        
        
        
        
        
    end
end

l2_analytic = -(8/3 + 56*y.^2 + 192*y.^4 + 160*y.^6);

figure
subplot(131)
surface(x,y,l2_a,'edgecolor','none')
title('X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')
colorbar
subplot(132)
surface(x,y,l2_b,'edgecolor','none')
title('X0^T*Q*X0-d.^2/m')
colorbar
subplot(133)
surface(x,y,l2_analytic,'edgecolor','none')
title('analytic')
colorbar


figure
subplot(121)
surface(x,y,abs((l2_a-l2_analytic))./abs(l2_analytic),'edgecolor','none')
title('rel error, X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')
colorbar()
subplot(122)
surface(x,y,abs((l2_b-l2_analytic))./abs(l2_analytic),'edgecolor','none')
title('rel error, X0^T*Q*X0-d.^2/m')
colorbar()



%{
s1_analytic = -(1+3*y.^2);
l1_analytic = 2 + 18*y.^2 + 24*y.^4;
l2_analytic = -(8/3 + 56*y.^2 + 192*y.^4 + 160*y.^6);


figure
subplot(131)
surface(x,y,s1_analytic,'edgecolor','none')
title('analytic s1')
colorbar()
subplot(132)
surface(x,y,l1_analytic,'edgecolor','none')
title('analytic lambda1')
colorbar()
subplot(133)
surface(x,y,l2_analytic,'edgecolor','none')
title('analytic lambda2')
colorbar()


figure
subplot(131)
surface(x,y,-s1_numerical(:,:,2),'edgecolor','none')
title('numerical s1')
colorbar()
subplot(132)
surface(x,y,l1_numerical(:,:,2),'edgecolor','none')
title('numerical lambda1')
colorbar()
subplot(133)
surface(x,y,l2_numerical(:,:,2),'edgecolor','none')
title('numerical lambda2')
colorbar()

figure
subplot(131)
surface(x,y,abs((s1_numerical(:,:,2))-s1_analytic)./abs(s1_analytic),'edgecolor','none')
title('rel error, s1')
colorbar()
subplot(132)
surface(x,y,abs((l1_numerical(:,:,2))-l1_analytic)./abs(l1_analytic),'edgecolor','none')
title('rel error, lambda1')
colorbar()
subplot(133)
surface(x,y,abs((l2_numerical(:,:,2))-l2_analytic)./abs(l2_analytic),'edgecolor','none')
title('rel error, lambda2')
colorbar()


k=0;
for T=0.001:0.001:0.8,
k=k+1;
T=-T;
sig_true=reshape(-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3)),[],1);
s1=reshape(-(1+3*y.^2),[],1);
lam1=reshape(2 + 18*y.^2 + 24*y.^4,[],1);
lam2=reshape(-(8/3 + 56*y.^2 + 192*y.^4 + 160*y.^6),[],1);

sig_0 = reshape(-s1,[],1);
sig_1 = reshape(-(-s1.^2 + 1/2*lam1) * T ,[],1);
sig_2 = reshape(-(4/3*s1.^3 - s1.*lam1 + 1/4*lam2) * T^2 ,[],1);

err0(k)=sqrt(mean((sig_true-sig_0).^2));
err1(k)=sqrt(mean((sig_true-(sig_0+sig_1)).^2));
err2(k)=sqrt(mean((sig_true-(sig_0+sig_1+sig_2)).^2));


sig_0_numerical = reshape(-s1_numerical(:,:,2),[],1);
sig_1_numerical = reshape(-(-s1_numerical(:,:,2).^2 + 1/2*l1_numerical(:,:,2)) * T ,[],1);
sig_2_numerical = reshape(-(4/3*s1_numerical(:,:,2).^3 - s1_numerical(:,:,2).*l1_numerical(:,:,2) + 1/4*l2_numerical(:,:,2)) * T^2 ,[],1);

err0_numerical(k)=sqrt(mean((sig_true-sig_0_numerical).^2));
err1_numerical(k)=sqrt(mean((sig_true-(sig_0_numerical+sig_1_numerical)).^2));
err2_numerical(k)=sqrt(mean((sig_true-(sig_0_numerical+sig_1_numerical+sig_2_numerical)).^2));



sig_0_plot(k) = mean(reshape(-s1,[],1));
sig_1_plot(k) = mean(reshape(-(-s1.^2 + 1/2*lam1) * T ,[],1));
sig_2_plot(k) = mean(reshape(-(4/3*s1.^3 - s1.*lam1 + 1/4*lam2) * T^2 ,[],1));
sig_0_numerical_plot(k) = mean(reshape(-s1_numerical(:,:,2),[],1));
sig_1_numerical_plot(k) = mean(reshape(-(-s1_numerical(:,:,2).^2 + 1/2*l1_numerical(:,:,2)) * T ,[],1));
sig_2_numerical_plot(k) = mean(reshape(-(4/3*s1_numerical(:,:,2).^3 - s1_numerical(:,:,2).*l1_numerical(:,:,2) + 1/4*l2_numerical(:,:,2)) * T^2 ,[],1));






%err0(k)=sqrt(immse(sig_true,sig_0));
%err1(k)=sqrt(immse(sig_true,sig_0+sig_1));
%err2(k)=sqrt(immse(sig_true,sig_0+sig_1+sig_2));
%rho0(k)=corr(      sig_true,sig_0);
%rho1(k)=corr(      sig_true,sig_0+sig_1);
tt(k)=T;
end
figure
hold on
plot(abs(tt),err0,'b')
plot(abs(tt),err1,'r')
plot(abs(tt),err2,'k')
plot(abs(tt),err0_numerical,'c--')
plot(abs(tt),err1_numerical,'m--')
plot(abs(tt),err2_numerical,'g--')
xlabel('|T|');
ylabel('Root mean-squared error');

figure
hold on
plot(abs(tt),sig_0_plot,'b')
plot(abs(tt),sig_1_plot,'r')
plot(abs(tt),sig_2_plot,'k')
plot(abs(tt),sig_0_numerical_plot,'c--')
plot(abs(tt),sig_1_numerical_plot,'m--')
plot(abs(tt),sig_2_numerical_plot,'g--')

%{
figure
hold on
plot(abs(t),err0_numerical)
plot(abs(t),err1_numerical,'m')
plot(abs(t),err2_numerical,'k')
xlabel('|T|');
ylabel('Root mean-squared error');


%}
%}