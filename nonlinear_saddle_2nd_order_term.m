close all
clear all
clc
%y=0:0.01:0.5;
%y=y(:);

xdim = 101;
ydim = 101;
xmin = 0
xmax = 0.5
x=linspace(xmin,xmax,xdim);
dx = x(2)-x(1);
y=linspace(xmin,xmax,ydim);
dy=y(2)-y(1);
[x,y]=meshgrid(x,y);

k=0;
for T=0.01:0.01:0.8,
k=k+1;
T=-T;
sig_0 = reshape((1+3*y.^2),[],1);
sig_true=reshape(-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3)),[],1);
err(k)=mean((sig_true-sig_0).^2);
%rho(k)=corr(sig_true,sig_0);
t(k)=T;
end

k=0;
for T=0.01:0.01:0.8,
k=k+1;
T=-T;
sig_0 = reshape((1+3*y.^2),[],1);
sig_1 = reshape(-3*T*y.^2.*(1+3*y.^2),[],1);
sig_true=reshape(-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3)),[],1);
err2(k)=mean((sig_true-(sig_0+sig_1)).^2);
%rho2(k)=corr(sig_true,sig_0+sig_1);
t(k)=T;
end
s1_analytic = 1+3*y.^2;
cor_analytic = 3*y.^2.*(1+y.^2);
figure
subplot(121)
surface(x,y,s1_analytic,'edgecolor','none')
title('analytic s1')
colorbar()
subplot(122)
surface(x,y,cor_analytic,'edgecolor','none')
title('analytic correction')
colorbar()

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

for t =1:3
    for i =1:ydim
        for j = 1:xdim
            Grad_v = [dudx(i,j,t), dudy(i,j,t);dvdx(i,j,t), dvdy(i,j,t)];
            Grad_D = [dDudx(i,j,t), dDudy(i,j,t); dDvdx(i,j,t), dDvdy(i,j,t)];
            S = 0.5*(Grad_v + Grad_v');
            B = 0.5*(Grad_D + Grad_D')+(Grad_v'*Grad_v);
            [V,D] = eig(S);
            if ~issorted(diag(D))
                [D,I] = sort(diag(D));
                V = V(:, I);
            end
            s1_numerical(i,j,t) = D(1,1);
            X1 = V(:,1);
            cor_numerical(i,j,t) = -s1_numerical(i,j,t).^2+0.5*(X1'*B*X1);
              
        end
    end
end

figure
subplot(121)
surface(x,y,-s1_numerical(:,:,2),'edgecolor','none')
title('numerical s1')
colorbar()
subplot(122)
surface(x,y,cor_numerical(:,:,2),'edgecolor','none')
title('numerical correction')
colorbar()


figure
subplot(121)
surface(x,y,abs((-s1_numerical(:,:,2))-s1_analytic),'edgecolor','none')
title('abs error, s1')
colorbar()
subplot(122)
surface(x,y,abs((cor_numerical(:,:,2))-cor_analytic),'edgecolor','none')
title('abs error, correction')
colorbar()

sig0_ap = reshape(-s1_numerical(:,:,2),[],1);
sig1_ap = reshape(cor_numerical(:,:,2),[],1);
k=0;
for T=0.0:0.01:0.8,
k=k+1;
T=-T;
sig_0 = reshape((1+3*y.^2),[],1);
sig_1 = reshape(-3*T*y.^2.*(1+3*y.^2),[],1);
sig_true=reshape(-(1/(2*T)*log(exp(4*T)./((1+y.^2)*exp(2*T) - y.^2).^3)),[],1);
err1(k)=mean((sig_true-sig_0).^2);
err2(k)=mean((sig_true-(sig_0+sig_1)).^2);
err1_num(k)=mean((sig_true-sig0_ap).^2);
err2_num(k)=mean((sig_true-(sig0_ap-T*sig1_ap)).^2);

%rho2(k)=corr(sig_true,sig_0+sig_1);
t(k)=T;
end
figure
hold on
plot(abs(t),err1,'m')
plot(abs(t),err2,'r')
plot(abs(t),err1_num,'c--')
plot(abs(t),err2_num,'b')