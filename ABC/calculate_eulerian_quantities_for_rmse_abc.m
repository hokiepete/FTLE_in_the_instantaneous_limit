close all
clear all
clc

dim = 96
tdim = 96
t0 = 0
tf = -0.84
x = linspace(0,2*pi,dim);
dx=x(2)-x(1);
dy=dx;
dz=dy;
twant = linspace(t0,tf,tdim);

[x,y,z]=meshgrid(x,x,x);

tdim = 3;
t0 = 0%0.5*pi
time = linspace(-0.01+t0,0.01+t0,tdim);
dt = time(2)-time(1);
for i =1:length(time)
    t=time(i);
    [uu,vv,ww]=abc_flow(x,y,z,t);
    u(:,:,:,i)=uu;
    v(:,:,:,i)=vv;
    w(:,:,:,i)=ww;
    
end
clear uu vv ww
[dudx,dudy,dudz,dudt] = gradient(u,dx,dy,dz,dt);
[dvdx,dvdy,dvdz,dvdt] = gradient(v,dx,dy,dz,dt);
[dwdx,dwdy,dwdz,dwdt] = gradient(w,dx,dy,dz,dt);

au = dudt+u.*dudx+v.*dudy+w.*dudz;
av = dvdt+u.*dvdx+v.*dvdy+w.*dvdz;
aw = dwdt+u.*dwdx+v.*dwdy+w.*dwdz;
clear dudt dvdt dwdt
[daudx,daudy,daudz,daudt] = gradient(au,dx,dy,dz,dt);
[davdx,davdy,davdz,davdt] = gradient(av,dx,dy,dz,dt);
[dawdx,dawdy,dawdz,dawdt] = gradient(aw,dx,dy,dz,dt);

ju = daudt+u.*daudx+v.*daudy+w.*daudz;
jv = davdt+u.*davdx+v.*davdy+w.*davdz;
jw = dawdt+u.*dawdx+v.*dawdy+w.*dawdz;
clear daudt davdt u v au av dawdt w aw
[djudx,djudy,djudz,djudt] = gradient(ju,dx,dy,dz,dt);
[djvdx,djvdy,djvdz,djvdt] = gradient(jv,dx,dy,dz,dt);
[djwdx,djwdy,djwdz,djwdt] = gradient(jw,dx,dy,dz,dt);
clear djudt djvdt djwdt
R = [0,-1;1,0];
for t =1:tdim
    t
    for i = 1:dim
        for j = 1:dim
            for k = 1:dim
                if ~isnan(ju(i,j,k,t))&&~isnan(jv(i,j,k,t))&&~isnan(jw(i,j,k,t))
                    Grad_v = [dudx(i,j,k,t), dudy(i,j,k,t), dudz(i,j,k,t);dvdx(i,j,k,t), dvdy(i,j,k,t), dvdz(i,j,k,t);dwdx(i,j,k,t), dwdy(i,j,k,t), dwdz(i,j,k,t)];
                    Grad_a = [daudx(i,j,k,t), daudy(i,j,k,t), daudz(i,j,k,t);davdx(i,j,k,t), davdy(i,j,k,t), davdz(i,j,k,t);dawdx(i,j,k,t), dawdy(i,j,k,t), dawdz(i,j,k,t)];
                    Grad_j = [djudx(i,j,k,t), djudy(i,j,k,t), djudz(i,j,k,t);djvdx(i,j,k,t), djvdy(i,j,k,t), djvdz(i,j,k,t);djwdx(i,j,k,t), djwdy(i,j,k,t), djwdz(i,j,k,t)];
                    S = 0.5*(Grad_v + Grad_v');
                    B = 0.5*(Grad_a + Grad_a')+(Grad_v'*Grad_v);
                    Q = 1./3.*(Grad_j + Grad_j')+(Grad_v'*Grad_a+Grad_a'*Grad_v);
                    [V,D] = eig(S);
                    if ~issorted(diag(D))
                        [D,I] = sort(diag(D));
                        V = V(:, I);
                    end
                    s1(i,j,k,t) = D(1,1);
                    X0 = V(:,1);
                    %xi(i,j,k,t,:) = V(:,1);
                    l1(i,j,k,t) = X0'*B*X0;
                    %
                    X1 = pinv((S-s1(i,j,k,t)*eye(size(B))))*(-((B-l1(i,j,k,t)*eye(size(B)))*X0));

                    l2(i,j,k,t) = X0'*Q*X0 + X0'*B*X1 - l1(i,j,k,t).*X0'*X1;
                    %db(i,j,k,t) = X0'*B*X1 - X0'*S*X1;
                    a1(i,j,k,t)=l2(i,j,k,t);
                    %Xp = R*X0;
                    %m = X0'*R'*(S-s1(i,j,k,t)*eye(size(S)))*R*X0;
                    %d = X0'*R'*B*X0;
                    %dd(i,j,k,t)=d;
                    %mm(i,j,k,t)=m;
                    %l2(i,j,k,t) = X0'*Q*X0-d.^2/m;%m/(X0'*Q*X0*m-d.^2);%X0'*Q*X0-d.^2/m;
                    %a2(i,j,k,t)=l2(i,j,k,t);
                    %X1 = V(:,1);
                    %l1(i,j,t) = X1'*B*X1;
                    %l2(i,j,t) = X1'*Q*X1;
                    %cor(i,j,t) = -s1(i,j,t).^2+0.5*(X1'*B*X1);

                else
                    s1(i,j,k,t) =nan;
                    l1(i,j,k,t) = nan;
                    l2(i,j,k,t) = nan;
                    %cor(i,j,t)=nan;
                end
            end
        end
    end
end

save abc_correction3rd s1 l1 l2
save abc_error_comparison s1 l1 a1 twant%a2% dd db
figure
vol3d_v2('cdata',a1)

%{
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
%

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

%}