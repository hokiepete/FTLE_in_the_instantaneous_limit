clear all
close all
clc
y=0:0.01:0.5;
dy = y(2)-y(1);
dx = dy;%y(2)-y(1);
dt=1;
ydim = length(y);
xdim=ydim;
[x,y]=meshgrid(y,y);
%y=y(:);


for t =1:3
    u(:,:,t) = x;
    v(:,:,t) = -y-y.^3;
end

[dudx,dudy,dudt] = gradient(u,dx,dy,1);
[dvdx,dvdy,dvdt] = gradient(v,dx,dy,1);

Du = dudt+u.*dudx+v.*dudy;
Dv = dvdt+u.*dvdx+v.*dvdy;

[dDudx,dDudy,dDudt] = gradient(Du,dx,dy,1);
[dDvdx,dDvdy,dDvdt] = gradient(Dv,dx,dy,1);

au = dudt;
av = dvdt;

[daudx,daudy,daudt] = gradient(au,dx,dy,dt);
[davdx,davdy,davdt] = gradient(av,dx,dy,dt);

[daudtdx,daudtdy,daudtdt] = gradient(daudt,dx,dy,dt);
[davdtdx,davdtdy,davdtdt] = gradient(davdt,dx,dy,dt);


Dau = daudt+u.*daudx+v.*daudy;
Dav = davdt+u.*davdx+v.*davdy;

[dDaudx,dDaudy,dDaudt] = gradient(Dau,dx,dy,dt);
[dDavdx,dDavdy,dDavdt] = gradient(Dav,dx,dy,dt);




for t =1:3
    for i =1:ydim
        for j = 1:xdim
            Grad_v = [dudx(i,j,t), dudy(i,j,t);dvdx(i,j,t), dvdy(i,j,t)];
            Grad_D = [dDudx(i,j,t), dDudy(i,j,t); dDvdx(i,j,t), dDvdy(i,j,t)];
            Grad_a = [daudx(i,j,t), daudy(i,j,t);davdx(i,j,t), davdy(i,j,t)];
            Grad_Da = [dDaudx(i,j,t), dDaudy(i,j,t); dDavdx(i,j,t), dDavdy(i,j,t)];
            Grad_dadt = [daudtdx(i,j,t), daudtdy(i,j,t); davdtdx(i,j,t), davdtdy(i,j,t)];
            S = 0.5*(Grad_v + Grad_v');
            B = 0.5*(Grad_D + Grad_D')+(Grad_v'*Grad_v);
            %Q = 0.5*(Grad_Da + Grad_Da')+(Grad_v'*Grad_a+Grad_a'*Grad_v);
            Q = 0.5*(Grad_dadt + Grad_dadt')+(Grad_v'*Grad_a+Grad_a'*Grad_v);
            [V,D] = eig(S);
            if ~issorted(diag(D))
                [D,I] = sort(diag(D));
                V = V(:, I);
            end
            s1_numerical(i,j,t) = D(1,1);
            X0 = V(:,1);
            l1_numerical(i,j,t) = X0'*B*X0;
            
            AA = (l1_numerical(i,j,t)*eye(size(B))-B)*X0;
            BB = S-s1_numerical(i,j,t)*eye(size(B));
            if BB(1,2)==0;
                X1(1,1)=AA(1);
            else
                error('uh oh b!=0')
            end
            
            if B(2,1)==0;
                X1(2,1)=AA(2);
            else
                error('uh oh c!=0')
            end
            %}
            
            %X1 = pinv(-((B-l1_numerical(i,j,t)*eye(size(B)))*X0)\(S-s1_numerical(i,j,t)*eye(size(B))));
            %X1=X1';
            %res(i,j,t,:)=-((B-l1_numerical(i,j,t)*eye(size(B)))*X0)-(S-s1_numerical(i,j,t)*eye(size(B)))*X1;
            l2_numerical(i,j,t) = X0'*Q*X0 + X0'*B*X1 - l1_numerical(i,j,t)*X0'*X1;
              
        end
    end
end


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