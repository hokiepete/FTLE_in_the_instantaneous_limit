
close all
clear all
clc

load high_res_vel
dx = xx(1,2)-xx(1,1)
dy = yy(2,1)-yy(1,1)
dt = time(2)-time(1)
[ydim,xdim,tdim]=size(uu);
%
[dudx,dudy,dudt] = gradient(uu,dx,dy,dt);
[dvdx,dvdy,dvdt] = gradient(vv,dx,dy,dt);

Du = dudt+uu.*dudx+vv.*dudy;
Dv = dvdt+uu.*dvdx+vv.*dvdy;

[dDudx,dDudy,dDudt] = gradient(Du,dx,dy,dt);
[dDvdx,dDvdy,dDvdt] = gradient(Dv,dx,dy,dt);

for t =1:tdim
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
                s1(i,j,t) = D(1,1);
                X1 = V(:,1);
                cor(i,j,t) = -s1_numerical(i,j,t).^2+0.5*(X1'*B*X1);

            else
                s1(i,j,t) =nan;
                cor(i,j,t)=nan;
            end
        end
    end
end

save correction_term s1 cor xx yy time