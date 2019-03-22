
close all
clear all
clc

load mv_vel_data
u(u==999)=nan;
v(v==999)=nan;

dt = 3600
dx = 200
dy = 200

[ydim,xdim,tdim] = size(u)
clear vv uu xx yy time
[dudx,dudy,dudt] = gradient(u,dx,dy,dt);
[dvdx,dvdy,dvdt] = gradient(v,dx,dy,dt);

au = dudt+u.*dudx+v.*dudy;
av = dvdt+u.*dvdx+v.*dvdy;
clear dudt dvdt
[daudx,daudy,daudt] = gradient(au,dx,dy,dt);
[davdx,davdy,davdt] = gradient(av,dx,dy,dt);

ju = daudt+u.*daudx+v.*daudy;
jv = davdt+u.*davdx+v.*davdy;
clear daudt davdt u v au av
[djudx,djudy,djudt] = gradient(ju,dx,dy,dt);
[djvdx,djvdy,djvdt] = gradient(jv,dx,dy,dt);
clear djudt djvdt
R = [0,-1;1,0];

for t =1:25%tdim
    t
    for i = 1:ydim
        for j = 1:xdim
            if ~isnan(ju(i,j,t))&&~isnan(jv(i,j,t))
                Grad_v = [dudx(i,j,t), dudy(i,j,t);dvdx(i,j,t), dvdy(i,j,t)];
                Grad_a = [daudx(i,j,t), daudy(i,j,t); davdx(i,j,t), davdy(i,j,t)];
                Grad_j = [djudx(i,j,t), djudy(i,j,t); djvdx(i,j,t), djvdy(i,j,t)];
                S = 0.5*(Grad_v + Grad_v');
                B = 0.5*(Grad_a + Grad_a')+(Grad_v'*Grad_v);
                Q = 0.5*(Grad_j + Grad_j')+(Grad_v'*Grad_a+Grad_a'*Grad_v);
                [V,D] = eig(S);
                if ~issorted(diag(D))
                    [D,I] = sort(diag(D));
                    V = V(:, I);
                end
                s1(i,j,t) = D(1,1);
                X0 = V(:,1);
                l1(i,j,t) = X0'*B*X0;
                X1 = -((B-l1(i,j,t)*eye(size(B)))*X0)\(S-s1(i,j,t)*eye(size(B)));
                X1=X1'/norm(X1);
                l2(i,j,t) = X0'*Q*X0 + X0'*B*X1 - X0'*S*X1;
                a1(i,j,t)=l2(i,j,t);
                %Xp = R*X0;
                m = X0'*R'*(S-s1(i,j,t)*eye(size(S)))*R*X0;
                d = X0'*R'*B*X0;
                l2(i,j,t) = X0'*Q*X0-d.^2/m;
                a2(i,j,t)=l2(i,j,t);
                %X1 = V(:,1);
                %l1(i,j,t) = X1'*B*X1;
                %l2(i,j,t) = X1'*Q*X1;
                %cor(i,j,t) = -s1(i,j,t).^2+0.5*(X1'*B*X1);
                a3(i,j,t) = norm(Q);
                a3(i,j,t) = norm(B);

            else
                s1(i,j,t) =nan;
                l1(i,j,t) = nan;
                l2(i,j,t) = nan;
                a1(i,j,t) = nan;
                a2(i,j,t) = nan;
                %cor(i,j,t)=nan;
            end
        end
    end
end

%save correction3rd s1 l1 l2
save error_comparison_terms s1 l1 a1 a2

figure
subplot(121)
surface(a1(:,:,end),'edgecolor','none')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')
colorbar
axis tight
subplot(122)
surface(a2(:,:,end),'edgecolor','none')
title('lambda2 = X0^T*Q*X0-d.^2/m')
axis tight
colorbar
