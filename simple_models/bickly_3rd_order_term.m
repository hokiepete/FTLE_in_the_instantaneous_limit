
close all
clear all
clc

ydim = 101;
xdim = ceil(2.5*ydim)
y=linspace(-4000,4000,ydim); %km
dy = y(2)-y(1);
x=linspace(0,20000,xdim); %km
dx = x(2)-x(1);
U0 = 62.66;%m s^-1
U0 = 24*3.6*U0 %km hr^-1
L = 1770;%km
A2 = 0.1;
A3 = 0.3;
c2 = 0.205*U0;
c3 = 0.461*U0;
re = 6371;%km
k = @(n) 2*n/re;
[x,y]=meshgrid(x,y);
%Fixed Reference Frame
tdim = 3;
t0 = 0%0.5*pi
time = linspace(-0.01+t0,0.01+t0,tdim);
dt = time(2)-time(1);
for i =1:length(time)
    t=time(i);
    dpdx = -k(3).*A3.*L.*U0.*((sech(y./L)).^2).*sin(k(3).*(x-c3.*t)) - k(2).*A2.*L.*U0.*((sech(y./L)).^2).*sin(k(2).*(x-c2.*t));
    dpdy = - U0.*((sech(y./L)).^2) - 2.*A3.*U0.*tanh(y./L).*((sech(y./L)).^2).*cos(k(3).*(x-c3.*t)) - 2.*A2.*U0.*tanh(y./L).*((sech(y./L)).^2).*cos(k(2).*(x-c2.*t));

    v(:,:,i) = dpdx;
    u(:,:,i) = -dpdy;

end

clear dpdx dpdy


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
b = zeros([ydim,xdim,tdim,2,2]);
for t =1:tdim
    t
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
                X1 = -(1./pinv((S-s1(i,j,t)*eye(size(B)))*(B-l1(i,j,t)*eye(size(B)))*X0));
                X1=X1';
                %{
                if sum(X1)~=0
                    X1=X1'/norm(X1);
                else
                    X1=X1';
                end
                %}
                l2(i,j,t) = X0'*Q*X0 + X0'*B*X1 - X0'*S*X1;
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

save bick_correction3rd s1 l1 l2
save bick_error_comparison s1 l1 a1 a2 dd db

figure
subplot(421)
surface(db(:,:,2),'edgecolor','none')
title('X0^T*B*X1 - X0^T*S*X1')
colorbar
axis tight
subplot(422)
surface(-dd(:,:,2).^2,'edgecolor','none')
title('-d^2')
axis tight
colorbar
subplot(423)
surface(mm(:,:,2),'edgecolor','none')
title('\mu')
axis tight
colorbar
subplot(424)
surface(-dd(:,:,2).^2./mm(:,:,2),'edgecolor','none')
title('-d^2/\mu')
axis tight
colorbar
subplot(425)
surface(-dd(:,:,2),'edgecolor','none')
title('-d')
axis tight
colorbar
subplot(426)
surface(-dd(:,:,2)./mm(:,:,2),'edgecolor','none')
title('-d/\mu')
axis tight
colorbar
subplot(427)
surface(-dd(:,:,2).^2./(mm(:,:,2)+1.5),'edgecolor','none')
title('-d^2/(\mu+1.5)')
axis tight
colorbar
subplot(428)
surface(-dd(:,:,2).^2./1.25,'edgecolor','none')
title('-d^2/(1.25)')
axis tight
colorbar

%
figure
subplot(121)
surface(a1(:,:,2),'edgecolor','none')
title('lambda2 = X0^T*Q*X0 + X0^T*B*X1 - X0^T*S*X1')
colorbar
axis tight
subplot(122)
surface(a2(:,:,2),'edgecolor','none')
title('lambda2 = X0^T*Q*X0-d.^2/m')
axis tight
colorbar


figure
subplot(421)
surface(x,y,-s1(:,:,2),'edgecolor','none')
title('-s1')
colorbar
subplot(422)
surface(x,y,-s1(:,:,2).^2+0.5*l1(:,:,2),'edgecolor','none')
title('correction')
colorbar
subplot(423)
surface(x,y,xi(:,:,2,1),'edgecolor','none')
title('xi\_x')
colorbar
subplot(424)
surface(x,y,xi(:,:,2,2),'edgecolor','none')
title('xi\_y')
colorbar
subplot(425)
surface(x,y,b(:,:,2,1,1),'edgecolor','none')
title('B\_xx')
colorbar
subplot(426)
surface(x,y,b(:,:,2,1,2),'edgecolor','none')
title('B\_xy')
colorbar
subplot(427)
surface(x,y,b(:,:,2,2,2),'edgecolor','none')
title('B\_yy')
colorbar
subplot(428)
surface(x,y,b(:,:,2,2,1),'edgecolor','none')
title('B\_yx')
colorbar
%}

figure
subplot(221)
surface(x,y,-s1(:,:,2),'edgecolor','none')
title('-s_1')
subplot(222)
surface(x,y,-s1(:,:,2).^2+0.5.*l1(:,:,2),'edgecolor','none')
title('correction')
subplot(223)
surface(x,y,xi(:,:,2,1),'edgecolor','none')
title('Xi_x')
subplot(224)
surface(x,y,xi(:,:,2,2),'edgecolor','none')
title('Xi_y')

figure
quiver(x(1:15:end,1:15:end),y(1:15:end,1:15:end),xi(1:15:end,1:15:end,2,1),xi(1:15:end,1:15:end,2,2))%,'edgecolor','none')

