
close all
clear all
clc

load high_res_vel
dt = time(2)-time(1)
dx = xx(1,2)-xx(1,1)
dy = yy(2,1)-yy(1,1)
u=uu; 
v=vv;
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
clear djudt djvdt djudt djvdt
R = [0,-1;1,0];
b = zeros([ydim,xdim,tdim,2,2]);
for t =1:3%tdim
    for i = 1:ydim
        for j = 1:xdim
            if ~isnan(ju(i,j,t))&&~isnan(jv(i,j,t))
                Grad_v = [dudx(i,j,t), dudy(i,j,t);dvdx(i,j,t), dvdy(i,j,t)];
                Grad_a = [daudx(i,j,t), daudy(i,j,t); davdx(i,j,t), davdy(i,j,t)];
                Grad_j = [djudx(i,j,t), djudy(i,j,t); djvdx(i,j,t), djvdy(i,j,t)];
                S = 0.5*(Grad_v + Grad_v');
                B = 0.5*(Grad_a + Grad_a')+(Grad_v'*Grad_v);
                b(i,j,t,:,:)=B;
                Q = 1./3.*(Grad_j + Grad_j')+(Grad_v'*Grad_a+Grad_a'*Grad_v);
                [V,D] = eig(S);
                if ~issorted(diag(D))
                    [D,I] = sort(diag(D));
                    V = V(:, I);
                end
                s1(i,j,t) = D(1,1);
                X0 = V(:,1);
                xi(i,j,t,:) = V(:,1);
                l1(i,j,t) = X0'*B*X0;
                %X1 = -((B-l1(i,j,t)*eye(size(B)))*X0)\(S-s1(i,j,t)*eye(size(B)));
                %X1 = -(1./pinv((S-s1(i,j,t)*eye(size(B))))*(B-l1(i,j,t)*eye(size(B)))*X0);
                X1 = pinv((S - s1(i,j,t)*eye(size(B))))*(-((B - l1(i,j,t)*eye(size(B)))*X0));
        
                %X1=X1';
                %{
                if sum(X1)~=0
                    X1=X1'/norm(X1);
                else
                    X1=X1';
                end
                %}
                l2(i,j,t) = X0'*Q*X0 + X0'*B*X1 - l1(i,j,t)*X0'*X1;
                db(i,j,t) = X0'*B*X1 - X0'*S*X1;
                a1(i,j,t)=l2(i,j,t);
                %Xp = R*X0;
                m = X0'*R'*(S-s1(i,j,t)*eye(size(S)))*R*X0;
                d = X0'*R'*B*X0;
                dd(i,j,t)=d;
                mm(i,j,t)=m;
                l2(i,j,t) = X0'*Q*X0-d.^2/m;
                a2(i,j,t)=l2(i,j,t);
                %X1 = V(:,1);
                %l1(i,j,t) = X1'*B*X1;
                %l2(i,j,t) = X1'*Q*X1;
                %cor(i,j,t) = -s1(i,j,t).^2+0.5*(X1'*B*X1);

            else
                s1(i,j,t) =nan;
                l1(i,j,t) = nan;
                l2(i,j,t) = nan;
                %cor(i,j,t)=nan;
            end
        end
    end
end
save wrf_correction3rd s1 l1 l2
save error_comparison s1 l1 a1 a2 dd db mm xi b

figure
subplot(321)
surface(db(:,:,2),'edgecolor','none')
title('X0^T*B*X1 - X0^T*S*X1')
colorbar
axis tight
subplot(322)
surface(-dd(:,:,2).^2,'edgecolor','none')
title('-d^2')
axis tight
colorbar
subplot(323)
surface(mm(:,:,2),'edgecolor','none')
title('\mu')
axis tight
colorbar
subplot(324)
surface(-dd(:,:,2).^2./mm(:,:,2),'edgecolor','none')
title('-d^2/\mu')
axis tight
colorbar
subplot(325)
surface(-dd(:,:,2),'edgecolor','none')
title('-d')
axis tight
colorbar
subplot(326)
surface(-dd(:,:,2)./mm(:,:,2),'edgecolor','none')
title('-d/\mu')
axis tight
colorbar


figure
subplot(121)
surface(a1(:,:,2),'edgecolor','none')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')
colorbar
axis tight
subplot(122)
surface(a2(:,:,2),'edgecolor','none')
title('lambda2 = X0^T*Q*X0-d.^2/m')
colorbar


figure
subplot(421)
surface(s1(:,:,2),'edgecolor','none')
title('s1')
axis tight
subplot(422)
surface(s1(:,:,2).^2-0.5*l1(:,:,2),'edgecolor','none')
title('correction')
axis tight
subplot(423)
surface(xi(:,:,2,1),'edgecolor','none')
title('xi\_x')
axis tight
subplot(424)
surface(xi(:,:,2,2),'edgecolor','none')
title('xi\_y')
axis tight
subplot(425)
surface(b(:,:,2,1,1),'edgecolor','none')
title('B\_xx')
axis tight
subplot(426)
surface(b(:,:,2,1,2),'edgecolor','none')
title('B\_xy')
axis tight
subplot(427)
surface(b(:,:,2,2,2),'edgecolor','none')
title('B\_yy')
axis tight
subplot(428)
surface(b(:,:,2,2,1),'edgecolor','none')
title('B\_yx')
axis tight
%}
