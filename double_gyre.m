close all
clear all
clc
%y=0:0.01:0.5;
%y=y(:);

xdim = 101;
ydim = 51;
xmin = 0
xmax = 2
ymin = 0;
ymax = 1;
x=linspace(xmin,xmax,xdim);
dx = x(2)-x(1);
y=linspace(ymin,ymax,ydim);
dy=y(2)-y(1);
[x,y]=meshgrid(x,y);

A = 0.1;
w = 0.2.*pi;
e = 0.25;

time = linspace(-0.01,0.01,3);
dt = time(2)-time(1);
for i =1:length(time)
    t=time(i);
    a = e.*sin(w.*t);
    b = 1-2.*e.*sin(w.*t);
    f = a.*x.^2+b.*x;
    dfdx = 2.*a.*x+b;    

    u(:,:,i) =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
    v(:,:,i) = pi.*A.*cos(pi.*f).*sin(y.*pi).*dfdx;
end

[dudx,dudy,dudt] = gradient(u,dx,dy,dt);
[dvdx,dvdy,dvdt] = gradient(v,dx,dy,dt);

Du = dudt+u.*dudx+v.*dudy;
Dv = dvdt+u.*dvdx+v.*dvdy;

[dDudx,dDudy,dDudt] = gradient(Du,dx,dy,dt);
[dDvdx,dDvdy,dDvdt] = gradient(Dv,dx,dy,dt);

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

%{
figure
surface(x,y,-s1_numerical(:,:,3)+1*cor_numerical(:,:,3),'edgecolor','none')
title('numerical')
colorbar()
%}
s1 = s1_numerical(:,:,2);
c1 = cor_numerical(:,:,2);
%load dg_ftle_data_long
load sigma_big

t0 = 0
tf = -0.96
time = linspace(t0,tf,101);

ftle=sigma;
ftle(1,:,:)=-s1;
n = length(time);
for i =1:n
    T=time(i);
    sig_approx = reshape(-s1-T*c1,[],1);
    sig_true = reshape(ftle(i,:,:),[],1);
    rmse_corrected(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2)
    denominator = den1*den2
    cor_corrected(i) = numerator./denominator;
    
    sig_approx = reshape(-s1,[],1);
    rmse_uncorrected(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2)
    denominator = den1*den2
    cor_uncorrected(i) = numerator./denominator;
    
    %cor(i) = corr(sig_approx,sig_true);
end

figure
subplot(121)
hold on
plot(-time,rmse_corrected,'b')
plot(-time,rmse_uncorrected,'r')
legend('-s1-T*corr','-s1')
ylabel('rmse')
xlabel('|T|')
subplot(122)
hold on
plot(-time,cor_corrected,'b')
plot(-time,cor_uncorrected,'r')
ylabel('correlation')
xlabel('|T|')
legend('-s1-T*corr','-s1')

save dg_plot_data rmse_corrected rmse_uncorrected time

%}