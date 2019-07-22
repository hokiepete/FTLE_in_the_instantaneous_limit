close all
clear all
clc

ncfile  = 'rtofs_glo_2ds_f000_1hrly_prog.nc'
ncid = netcdf.open(ncfile);
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)

y = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Y')));
y = y(1:end-1);
x = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'X')));

lat = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Latitude')));
lat = lat(:,1:end-1);
lon = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Longitude')));
lon = lon(:,1:end-1);
u = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u_velocity')));
u = u(:,1:end-1);
u(u>1e10)=nan;
%u(lon>365)=nan;
v = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v_velocity')));
v = v(:,1:end-1);
v(v>1e10)=nan;
%v(lon>365)=nan;
netcdf.close(ncid);
%lon = mod(lon-min(min(lon)),365);
[x,y]=ndgrid(x,y);
figure
surface(x,y,sqrt(u.^2+v.^2),'edgecolor','none');
colorbar
axis tight
title('speed')
saveas(gcf,'speed_ocean.fig')
figure
surface(lon,lat,u,'edgecolor','none');
colorbar
axis tight
title('u')
saveas(gcf,'u_ocean.fig')
figure
surface(lon,lat,v,'edgecolor','none');
colorbar
axis tight
title('v')
saveas(gcf,'v_ocean.fig')


%{
u(abs(lat)>66) = [];
%u = reshape(u,4500,[]);
v(abs(lat)>66) = [];
%v = reshape(v,4500,[]);
lon(abs(lat)>66) = [];
%lon = reshape(lon,4500,[]);
lat(abs(lat)>66) = [];
%lat = reshape(lat,4500,[]);
%}

%{
saveas(gcf,'ocean_speed.png')
close gcf
figure('visible','off');
plot(reshape(x,[],1),reshape(y,[],1),'.')
saveas(gcf,'scatter_x_y.png')
close gcf
figure('visible','off');
plot(reshape(lon,[],1),reshape(lat,[],1),'.')
saveas(gcf,'scatter_lon_lat.png')
%}