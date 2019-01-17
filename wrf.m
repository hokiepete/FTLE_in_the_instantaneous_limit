%{
close all
clear all
clc
%
dt=1 %hr
dx=12 %kms
dy=12 %kms


ncfile='ftle_80m.nc';
ncid=netcdf.open(ncfile,'NC_NOWRITE');

%{
Grab NC data
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
%
%Loop through NC data infomation
for i = 1:nvars
    i-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i-1)
end
%}

time_ftle = netcdf.getVar(ncid,0); %days
time_ftle=24*time_ftle(1:6:end); %hrs
lon = netcdf.getVar(ncid,1,'double');
lat = netcdf.getVar(ncid,2,'double');
ftle = netcdf.getVar(ncid,5); %days^{-1}
ftle=permute(squeeze(ftle(1,:,:,1:6:end)),[2,1,3]);
ftle(ftle==999999)=nan;
ftle=1/24*ftle; %hrs^{-1}
%
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

for t =1:tdim
    t
    for i =1:ydim
        for j = 1:xdim
            if ~isnan(Du(i,j,t))|~isnan(Dv(i,j,t))
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
%ftle(s1_numerical==nan)=nan;
%ftle(cor_numerical==nan)=nan;
%
figure
subplot(121)
surface(-s1_numerical(:,:,3),'edgecolor','none')
title('numerical s1')
colorbar()
subplot(122)
surface(cor_numerical(:,:,3),'edgecolor','none')
title('numerical correction')
colorbar()
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
time = time_ftle;
n = length(time);
for i =1:n
    T=-time(i);
    ftle_t = squeeze(ftle(:,:,end+1-i))
    sig_true = reshape(ftle_t,[],1);
    sig_approx = reshape(-s1-T*c1,[],1);
    ind = ~isnan(sig_true) & ~isnan(sig_approx) ;
    sig_true = sig_true(ind);
    sig_approx = sig_approx(ind);
    rmse_corrected(i) = sqrt(mean((sig_approx-sig_true).^2));
    sa_bar = mean(sig_approx);
    st_bar = mean(sig_true);
    numerator = sum(sig_approx.*sig_true)-(n*sa_bar*st_bar);
    den1 = sqrt(sum(sig_approx.^2)-n*sa_bar.^2);
    den2 = sqrt(sum(sig_true.^2)-n*st_bar.^2)
    denominator = den1*den2
    cor_corrected(i) = numerator./denominator;
    
    sig_approx = reshape(-s1,[],1);
    sig_approx = sig_approx(ind);
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
plot(time,rmse_corrected,'b')
plot(time,rmse_uncorrected,'r')
legend('-s1-T*corr','-s1')
ylabel('rmse')
xlabel('|T|')
subplot(122)
hold on
plot(time,cor_corrected,'b')
plot(time,cor_uncorrected,'r')
ylabel('correlation')
xlabel('|T|')
legend('-s1-T*corr','-s1')

%}