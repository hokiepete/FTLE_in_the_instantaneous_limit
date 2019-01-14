clear all
close all
clc
len = 96
t0 = 0
x = linspace(0,2*pi,len);
dx=x(2)-x(1);
[x,y,z] = meshgrid(x,x,x);
%
[u,v,w] = abc_flow(x,y,z,t0);
[dudx,dudy,dudz] = gradient(u,dx);
[dvdx,dvdy,dvdz] = gradient(v,dx);
[dwdx,dwdy,dwdz] = gradient(w,dx);

for i = 1:len
    for j = 1:len
        for k = 1:len
            gradU = [dudx(i,j,k),dudy(i,j,k),dudz(i,j,k);
                     dvdx(i,j,k),dvdy(i,j,k),dvdz(i,j,k);
                     dwdx(i,j,k),dwdy(i,j,k),dwdz(i,j,k)];
            S = 0.5*(gradU+gradU');
            s1(i,j,k) = min(eig(S));
        end
    end
end
save('s1_small.mat','s1');
figure
vol3d_v2('cdata',s1)
xlabel('x')
ylabel('y')
zlabel('z')
colorbar

%quiver3(x,y,z,u,v,w)
%vol3d_v2('cdata',u)
%}
%{

x = linspace(0,2*pi,len);
dx=x(2)-x(1);
twant = fliplr(linspace(t0,tf,len));

for i = 1:len
    i
    for j = 1:len
        for k = 1:len
            y0=[x(i,j,k),y(i,j,k),z(i,j,k)];
            [t,yout] = ode45(@abc_int,[tf,t0],y0);
            fx(:,i,j,k) = interp1(t,yout(:,1),twant);
            fy(:,i,j,k) = interp1(t,yout(:,2),twant);
            fz(:,i,j,k) = interp1(t,yout(:,3),twant);
        end
    end
end
sigma = zeros([len,len,len,len]);
for t = 1:len
    t
    %backward-time integration starts at tf and ends at t0
    %need to skip tf to avoid division by zero
    %set FTLE to zero instead
    if twant(t) == tf
        continue
    else
        [dfxdx,dfxdy,dfxdz] = gradient(squeeze(fx(t,:,:,:)),dx);
        [dfydx,dfydy,dfydz] = gradient(squeeze(fy(t,:,:,:)),dx);
        [dfzdx,dfzdy,dfzdz] = gradient(squeeze(fz(t,:,:,:)),dx);
        for i = 1:len
            for j = 1:len
                for k = 1:len
                    gradF = [dfxdx(i,j,k),dfxdy(i,j,k),dfxdz(i,j,k);
                             dfydx(i,j,k),dfydy(i,j,k),dfydz(i,j,k);
                             dfzdx(i,j,k),dfzdy(i,j,k),dfzdz(i,j,k)];
                    C = gradF'*gradF;
                    lambda = max(eig(C));
                    sigma(t,i,j,k) = 1/(2*abs(twant(t)-tf))*log(lambda);


                end
            end
        end
    end
end

save('output.mat', 'sigma');
%}
%{
figure
vol3d_v2('cdata',squeeze(sigma(1,:,:,:)))
xlabel('x')
ylabel('y')
zlabel('z')
colorbar

figure
vol3d_v2('cdata',squeeze(sigma(end,:,:,:)))
xlabel('x')
ylabel('y')
zlabel('z')
colorbar

%}


clear all
close all
clc
len = 96
tlen = 9
t0 = 0
tf = -1
x = linspace(0,2*pi,len);
dx=x(2)-x(1);
twant = linspace(t0,tf,tlen);
[x,y,z] = meshgrid(x,x,x);
for i = 1:len
    i
    for j = 1:len
        for k = 1:len
            y0=[x(i,j,k),y(i,j,k),z(i,j,k)];
            [t,yout] = ode45(@abc_int,[t0,tf],y0);
            fx(:,i,j,k) = interp1(t,yout(:,1),twant);
            fy(:,i,j,k) = interp1(t,yout(:,2),twant);
            fz(:,i,j,k) = interp1(t,yout(:,3),twant);
        end
    end
end
sigma = zeros([tlen,len,len,len]);
for t = 1:tlen
    t
    %need to skip t0 to avoid division by zero
    %set FTLE to zero instead
    if twant(t) == t0
        continue
    else
        [dfxdx,dfxdy,dfxdz] = gradient(squeeze(fx(t,:,:,:)),dx);
        [dfydx,dfydy,dfydz] = gradient(squeeze(fy(t,:,:,:)),dx);
        [dfzdx,dfzdy,dfzdz] = gradient(squeeze(fz(t,:,:,:)),dx);
        for i = 1:len
            for j = 1:len
                for k = 1:len
                    gradF = [dfxdx(i,j,k),dfxdy(i,j,k),dfxdz(i,j,k);
                             dfydx(i,j,k),dfydy(i,j,k),dfydz(i,j,k);
                             dfzdx(i,j,k),dfzdy(i,j,k),dfzdz(i,j,k)];
                    C = gradF'*gradF;
                    lambda = max(eig(C));
                    sigma(t,i,j,k) = 1/(2*abs(twant(t)-t0))*log(lambda);


                end
            end
        end
    end
end

save('sigma_small.mat', 'sigma');
