
close all
clear all
clc

dt=1 %hr
%dx=3.4 %kms
%dy=3.4 %kms

%nums = linspace(3,4,101);
dx = 3.483;
dy=dx
%{
dt=3600 %s
dx=12*1000 %m
dy=12*1000 %m
%}
ncfile='hosiendata_wind_velocity.nc';
ncid=netcdf.open(ncfile,'NC_NOWRITE');
time_U = netcdf.getVar(ncid,0);
u = netcdf.getVar(ncid,3); %m/s
u = permute(u,[2,1,3]);
u(u==999)=nan;
u=3.6*u; %km/hr
v = netcdf.getVar(ncid,4); %m/s
v = permute(v,[2,1,3]);
v(v==999)=nan;
v=3.6*v; %km/hr
[ydim,xdim,tdim]=size(u);
%
[dudx,dudy,dudt] = gradient(u,dx,dy,dt);
[dvdx,dvdy,dvdt] = gradient(v,dx,dy,dt);

Du = dudt+u.*dudx+v.*dudy;
Dv = dvdt+u.*dvdx+v.*dvdy;

[dDudx,dDudy,dDudt] = gradient(Du,dx,dy,dt);
[dDvdx,dDvdy,dDvdt] = gradient(Dv,dx,dy,dt);

for t =25:length(time_U)
    t
    for i =1:ydim
        for j = 1:xdim
            if ~isnan(Du(i,j,t))&&~isnan(Dv(i,j,t))
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

            else
                s1_numerical(i,j,t) =nan;
                cor_numerical(i,j,t)=nan;
            end
        end
    end
end

ncfile='ftle_80m_small.nc';
ncid=netcdf.open(ncfile,'NC_NOWRITE');
%{
%Grab NC data
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
%
%Loop through NC data infomation
for i = 1:nvars
    i-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i-1)
end
%}

time_ftle = netcdf.getVar(ncid,0); %days
%time_ftle=24*time_ftle(1:6:end); %hrs
time_ftle=24*time_ftle; %hrs
lon = netcdf.getVar(ncid,1,'double');
lat = netcdf.getVar(ncid,2,'double');
ftle = netcdf.getVar(ncid,4); %days^{-1}
%ftle=permute(squeeze(ftle(1,:,:,1:6:end)),[2,1,3]);
ftle=permute(squeeze(ftle(1,:,:,:)),[2,1,3]);
ftle(ftle==999999)=nan;
%ftle=1/(24*3.4)*ftle; %hrs^{-1}
ftle=1/(24)*ftle; %hrs^{-1}
[ydim,xdim,tdim]=size(ftle)
%
%ftle(s1_numerical==nan)=nan;
%ftle(cor_numerical==nan)=nan;
%}

%
%{
figure
surface(x,y,-s1_numerical(:,:,3)+1*cor_numerical(:,:,3),'edgecolor','none')
title('numerical')
colorbar()
%}
%}
s1 = s1_numerical(:,:,end);
c1 = cor_numerical(:,:,end);
ftle(:,:,end)=-s1;
time = 22-time_ftle;
n = length(time);
for i =1:n
    T=time(i);
    ftle_t = squeeze(ftle(:,:,end+1-i));
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
close all
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
size(time);
size(rmse_corrected);
%subplot(121)
%
t_want=linspace(0,0.7,101);
int = 1
yy = interp1(-time(1:int:end),rmse_corrected(1:int:end),t_want,'spline')
xx = interp1(-time(1:int:end),rmse_uncorrected(1:int:end),t_want,'spline')

%yy = interp1(-time,rmse_corrected,t_want,'spline')
fig=figure
hold on
plot(-time(1:end),rmse_corrected(1:end),'b.-')
plot(t_want,yy,'b.-')
plot(-time(1:end),rmse_uncorrected(1:end),'r-')
plot(t_want,xx,'r.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE hr^{-1}')
xlabel('|T| hr')
xlim([0,0.7])
%    
save wrf_plot_data rmse_corrected rmse_uncorrected time
%{
subplot(122)
hold on
plot(time,cor_corrected,'b')
plot(time,cor_uncorrected,'r')
ylabel('correlation')
xlabel('|T|')
legend('-s1-T*corr','-s1','Location','southwest')
%saveas(fig,sprintf('%1.1f.fig',Q))
%
weight(j) = Q
cor_min(j) = min(rmse_corrected)
uncor_min(j) = min(rmse_uncorrected)
end
figure;hold on;plot(weight,cor_min,'b-');plot(weight,uncor_min,'r-')
figure
subplot(121)
surface(-s1,'edgecolor','none')
colorbar
subplot(122)
surface(ftle(:,:,1),'edgecolor','none')
colorbar

%}