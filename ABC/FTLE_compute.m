
clear all
close all
clc

len = 96;
x = linspace(0,2*pi,len);
dx = x(2)-x(1);
[x,y,z]=meshgrid(x,x,x);

n=25;
tend=-1
t = linspace(0,tend,n);
dt = t(2)-t(1);
for i = 1:len
    i
    for j=1:len
        for k=1:len
            y0=[x(i,j,k),y(i,j,k),z(i,j,k)];
            [tout,yout] = ode45(@abc_int,[0,tend],y0);
            fx(i,j,k,:) = interp1(tout,yout(:,1),t);
            fy(i,j,k,:) = interp1(tout,yout(:,2),t);
            fz(i,j,k,:) = interp1(tout,yout(:,3),t);

        end
    end
end

[dfxdx,dfxdy,dfxdz,dfxdt] = gradient(fx,dx,dx,dx,dt);
[dfydx,dfydy,dfydz,dfydt] = gradient(fy,dx,dx,dx,dt);
[dfzdx,dfzdy,dfzdz,dfzdt] = gradient(fz,dx,dx,dx,dt);

clear fx fy fz dfxdt dfydt dfzdt

sigma=zeros(len,len,len,n);
for tt = 2:n
    for i = 1:len
    i
        for j=1:len
            for k=1:len
                grad_F = [dfxdx(i,j,k,tt),dfxdy(i,j,k,tt),dfxdz(i,j,k,tt);dfydx(i,j,k,tt),dfydy(i,j,k,tt),dfydz(i,j,k,tt);dfzdx(i,j,k,tt),dfzdy(i,j,k,tt),dfzdz(i,j,k,tt)];
                C = grad_F'*grad_F;
                lambda = max(eig(C));
                sigma(i,j,k,tt) = 1/(2*abs(t(tt)-0))*log(lambda);
                
            end
        end
    end
end

save abc_ftle sigma t

vol3d_v2('cdata',sigma(:,:,:,end))