clear all
close all
clc
%{
tlen = 101
t0 = 0
tf = -0.96;%-2.5
time = linspace(t0,tf,tlen);
%
%time(2)-time(1);
x = linspace(0,2,101);
dx = x(2)-x(1);
y = linspace(0,1,51);
dy = y(2)-y(1);
dz=dy;
[x,y,z]=meshgrid(x,y,y);
gradient_time=linspace(-0.01,0.01,3);
for i =1:3
    [u(:,:,:,i),v(:,:,:,i),w(:,:,:,i)]=dg_flow(x,y,z,gradient_time(i));
end
dt = gradient_time(2)-gradient_time(1);
clear x y z gradient_time

[ydim,xdim,zdim,tdim]=size(u);

[dudx,dudy,dudz,dudt] = gradient(u,dx,dy,dz,dt);
[dvdx,dvdy,dvdz,dvdt] = gradient(v,dx,dy,dz,dt);
[dwdx,dwdy,dwdz,dwdt] = gradient(w,dx,dy,dz,dt);

Du = dudt+u.*dudx+v.*dudy+w.*dudz;
Dv = dvdt+u.*dvdx+v.*dvdy+w.*dvdz;
Dw = dwdt+u.*dwdx+v.*dwdy+w.*dwdz;
clear u v w dudt dvdt dwdt

[dDudx,dDudy,dDudz,dDudt] = gradient(Du,dx,dy,dz,dt);
[dDvdx,dDvdy,dDvdz,dDvdt] = gradient(Dv,dx,dy,dz,dt);
[dDwdx,dDwdy,dDwdz,dDwdt] = gradient(Dw,dx,dy,dz,dt);
clear Du Dv Dw
for t =1:3%length(time)
    t
    for i =1:ydim
        for j = 1:xdim
            for k = 1:zdim
                %if ~isnan(Du(i,j,k,t))&&~isnan(Dv(i,j,k,t))&&~isnan(Dw(i,j,k,t))
                    Grad_v = [dudx(i,j,k,t), dudy(i,j,k,t), dudz(i,j,k,t);dvdx(i,j,k,t), dvdy(i,j,k,t), dvdz(i,j,k,t); dwdx(i,j,k,t), dwdy(i,j,k,t), dwdz(i,j,k,t)];
                    Grad_D = [dDudx(i,j,k,t), dDudy(i,j,k,t), dDudz(i,j,k,t);dDvdx(i,j,k,t), dDvdy(i,j,k,t), dDvdz(i,j,k,t); dDwdx(i,j,k,t), dDwdy(i,j,k,t), dDwdz(i,j,k,t)];
                    S = 0.5*(Grad_v + Grad_v');
                    B = 0.5*(Grad_D + Grad_D')+(Grad_v'*Grad_v);
                    [V,D] = eig(S);
                    if ~issorted(diag(D))
                        [D,I] = sort(diag(D));
                        V = V(:, I);
                    end
                    s1(i,j,k,t) = D(1,1);
                    X1 = V(:,1);
                    correction(i,j,k,t) = -s1(i,j,k,t).^2+0.5*(X1'*B*X1);

                %else
                    %s1_numerical(i,j,k,t) =nan;
                    %cor_numerical(i,j,k,t)=nan;
                %end
            end
        end
    end
end
s1 = s1(:,:,:,2);
correction = correction(:,:,:,2);
clear dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz dDudx dDudy dDudz dDvdx dDvdy dDvdz dDwdx dDwdy dDwdz
save eulerian_rmse_data s1 correction% time
%}
load sigma_big
load eulerian_rmse_data
tlen = 101
t0 = 0
tf = -0.96;%-2.5
time = linspace(t0,tf,tlen);

c1 = correction;
ftle = sigma;
ftle(1,:,:,:)=-s1;
n = length(time);
for i =1:n
    T=time(i);
    ftle_t = squeeze(ftle(i,:,:,:));
    sig_true = reshape(ftle_t,[],1);
    sig_approx = reshape(-s1-T*c1,[],1);
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
    cor_corrected(i) = numerator./denominator;

    sig_approx = reshape(-s1,[],1);
    sig_approx = sig_approx(ind);
    rmse_uncorrected(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2);
    denominator = den1*den2;
    cor_uncorrected(i) = numerator./denominator;

    %cor(i) = corr(sig_approx,sig_true);
end
%time=time(2:27)
start = 1;
stop = 101;
close all
%time = time*60;
time=-time(start:stop);
rmse_corrected=rmse_corrected(start:stop);
rmse_uncorrected=rmse_uncorrected(start:stop);
cor_corrected=cor_corrected(start:stop);
cor_uncorrected=cor_uncorrected(start:stop);
size(time);
size(rmse_corrected);
%subplot(121)
hold on
plot(time,rmse_corrected,'b.-')
plot(time,rmse_uncorrected,'r.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE hr^{-1}')
xlabel('|T| hr')
xlim([0,0.7])
save dg3d_plot_data rmse_corrected rmse_uncorrected time
