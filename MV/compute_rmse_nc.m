close all
clear all
clc

load error_comparison_terms

ncfile='mv_ftle.nc';
ncfile='windage_ftle.nc';
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

time_ftle = netcdf.getVar(ncid,0)-74; %days
time_ftle=24*3600*time_ftle; %s
lon = netcdf.getVar(ncid,1,'double');
lat = netcdf.getVar(ncid,2,'double');
sigma = netcdf.getVar(ncid,5); %days^{-1}
sigma=permute(squeeze(sigma(1,:,:,:)),[2,1,3]);
sigma(sigma==999999)=nan;
sigma=1/(24*3600)*sigma; %s^{-1}
[ydim,xdim,tdim]=size(sigma)
%{
figure
surface(ftle(:,:,1),'edgecolor','none')
%}

s1 = s1(:,:,end);
l1 = l1(:,:,end);
a1 = a1(:,:,end);
a2 = a2(:,:,end);

%{
i=101;
s1 = s1((i+1):end-i,(i+1):end-i,1);
l1 = l1((i+1):end-i,(i+1):end-i,1);
a1 = a1((i+1):end-i,(i+1):end-i,1);
a2 = a2((i+1):end-i,(i+1):end-i,1);
sigma = sigma((i+1):end-i,(i+1):end-i,:);
%}

sigma(:,:,end)=-s1;
n = length(time_ftle);
for i =1:n
    T = time_ftle(i);
    ftle_t = squeeze(sigma(:,:,end-i+1));
    sig_true = reshape(ftle_t,[],1);
    sig_approx = reshape(-s1-T*(-s1.^2+0.5*l1+T*(4/3*s1.^3-s1.*l1+0.25*a1)),[],1);
    sig_approx2 = reshape(-s1-T*(-s1.^2+0.5*l1+T*(4/3*s1.^3-s1.*l1+0.25*a2)),[],1);
    
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
    
    sig_approx = reshape(-s1-T*(-s1.^2+0.5*l1),[],1);
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
plot(time_ftle,rmse1,'r.-')
plot(time_ftle,rmse2,'b.-')
%plot(time_ftle,rmsea1,'k.-')
legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')

subplot(122)
hold on
plot(time_ftle,rmse1,'r.-')
plot(time_ftle,rmse2,'b.-')
%plot(time_ftle,rmsea2,'k.-')
legend('-s1','-s1-O(T)','-s1-O(T^2)','Location','northwest')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
title('lambda2 = X0^T*Q*X0-d.^2/m')

%{
figure
subplot(121)
hold on
plot(time_ftle,rmse1,'r.-')
plot(time_ftle,rmse2,'b.-')
plot(time_ftle,rmsea1,'k.-')
legend('-s1-T*corr','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
subplot(122)
hold on
plot(time_ftle,rmse1,'r.-')
plot(time_ftle,rmse2,'b.-')
plot(time_ftle,rmsea2,'k.-')
legend('-s1-T*cor','-s1','Location','southeast')
ylabel('RMSE s^{-1}')
xlabel('|T| s')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis('equal')
%}
rmse3 = rmsea1;
time = T;
save wrf_plot_data rmse1 rmse2 rmse3 time


