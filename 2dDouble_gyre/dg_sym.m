close all
clear all
clc

%{

%%% Arbitrary f(x,t)
syms u(x,y,t) v(x,y,t) f(x,t) e w A

assume(A,'real')
assume(e,'real')
assume(w,'real')
assume(x,'real')
assume(y,'real')
assume(t,'real')

u(x,y,t) =-pi*A*sin(pi*f)*cos(y*pi);    
v(x,y,t) = pi*A*cos(pi*f)*sin(y*pi)*diff(f,x);

dudt = diff(u,t);
dudx = diff(u,x);
dudy = diff(u,y);

dvdt = diff(v,t);
dvdx = diff(v,x);
dvdy = diff(v,y);

au = dudt + u*dudx + v*dudy;
av = dvdt + u*dvdx + v*dvdy;

daudt = diff(au,t);
daudx = diff(au,x);
daudy = diff(au,y);

davdt = diff(av,t);
davdx = diff(av,x);
davdy = diff(av,y);

ju = daudt+u.*daudx+v.*daudy;
jv = davdt+u.*davdx+v.*davdy;

djudt = diff(ju,t);
djudx = diff(ju,x);
djudy = diff(ju,y);

djvdt = diff(jv,t);
djvdx = diff(jv,x);
djvdy = diff(jv,y);


grad_v = [dudx,dudy;dvdx,dvdy];
grad_a = [daudx,daudy;davdx,davdy];
grad_j = [djudx,djudy;djvdx,djvdy];

S = 0.5*(grad_v+grad_v');
B = 0.5*(grad_a+grad_a')+grad_v'*grad_v;
Q = 0.5*(grad_j + grad_j')+(grad_v'*grad_a+grad_a'*grad_v);

%}

%

%%% Known f(x,t)

syms u(x,y,t) v(x,y,t) f(x,t)% e w A
%assume(A,'real')
%assume(e,'real')
%assume(w,'real')
assume(x,'real')
assume(y,'real')
assume(t,'real')

A=0.1;
w=0.2*pi;
e=0.25;

f(x,t) = (e.*sin(w.*t)).*x.^2+(1-2.*e.*sin(w.*t)).*x;
u(x,y,t) =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
v(x,y,t) = pi.*A.*cos(pi.*f).*sin(y.*pi).*diff(f,x);
dudt = diff(u,t);
dudx = diff(u,x);
dudy = diff(u,y);

dvdt = diff(v,t);
dvdx = diff(v,x);
dvdy = diff(v,y);

au = dudt + u*dudx + v*dudy;
av = dvdt + u*dvdx + v*dvdy;

daudt = diff(au,t);
daudx = diff(au,x);
daudy = diff(au,y);

davdt = diff(av,t);
davdx = diff(av,x);
davdy = diff(av,y);

ju = daudt+u.*daudx+v.*daudy;
jv = davdt+u.*davdx+v.*davdy;

djudt = diff(ju,t);
djudx = diff(ju,x);
djudy = diff(ju,y);

djvdt = diff(jv,t);
djvdx = diff(jv,x);
djvdy = diff(jv,y);

grad_v = [dudx,dudy;dvdx,dvdy];
grad_a = [daudx,daudy;davdx,davdy];
grad_j = [djudx,djudy;djvdx,djvdy];

%{
grad_v = [dudx,dvdx;dudy,dvdy];
grad_a = [daudx,davdx;daudy,davdy];
grad_j = [djudx,djvdx;djudy,djvdy];
%}

S = 0.5*(grad_v+grad_v');
B = 0.5*(grad_a+grad_a')+grad_v'*grad_v;
Q = 0.5*(grad_j + grad_j')+(grad_v'*grad_a+grad_a'*grad_v);

xdim = 301
ydim =151
xx = linspace(0,2,xdim);
yy = linspace(0,1,ydim);
[xx,yy]=meshgrid(xx,yy);
R = [0,-1;1,0];
tt = 0;
for i = 1:ydim
    for j = 1:xdim
        i,j
        SS = reshape(double(subs(S,[x,y,t],[xx(i,j),yy(i,j),tt])),2,2);
        BB = reshape(double(subs(B,[x,y,t],[xx(i,j),yy(i,j),tt])),2,2);
        QQ = reshape(double(subs(Q,[x,y,t],[xx(i,j),yy(i,j),tt])),2,2);
        [V D] = eig(SS);
        if ~issorted(diag(D))
            [D,I] = sort(diag(D));
            V = V(:, I);
        end
        lambda_0(i,j) = D(1,1);
        X0 = V(:,1);
        lambda_1(i,j) = X0'*BB*X0;
        
        %First lambda_2 method calculating Xi_1
        X1 = -((BB - lambda_1(i,j)*eye(size(BB)))*X0) \ (SS - lambda_0(i,j)*eye(size(SS)));
        if sum(X1)~=0
            X1=X1'/norm(X1);
        else
            X1=X1';
        end
        %lambda_2_first(i,j) = X0'*QQ*X0 + X0'*BB*X1 - X0'*SS*X1;
        lambda_2_first(i,j) = X0'*QQ*X0 + X0'*BB*X1 - lambda_1(i,j).*X0'*X1;
        
        %Second lambda_2 method bypassing Xi_1
        mu = X0'*R'*(SS-lambda_0(i,j)*eye(size(SS)))*R*X0;
        d = X0'*R'*BB*X0;
        if mu~=0
            lambda_2_second(i,j) = X0'*QQ*X0 - d.^2./mu;
        else
            lambda_2_second(i,j)=nan;
        end
    end
end
figure
subplot(221)
surface(xx,yy,lambda_0,'edgecolor','none')
colorbar
subplot(222)
surface(xx,yy,lambda_1,'edgecolor','none')
colorbar
subplot(223)
surface(xx,yy,lambda_2_first,'edgecolor','none')
colorbar
subplot(224)
surface(xx,yy,lambda_2_second,'edgecolor','none')
colorbar

save analytic_lambda_terms xx yy lambda_0 lambda_1 lambda_2_first lambda_2_second


figure
subplot(221)
surface(xx,yy,check,'edgecolor','none')
colorbar
subplot(222)
surface(xx,yy,heck,'edgecolor','none')
colorbar
subplot(223)
surface(xx,yy,eck,'edgecolor','none')
colorbar
subplot(224)
surface(xx,yy,ck,'edgecolor','none')
colorbar
%}

%{

%%% derivatives of f(x,t)

syms f(x,t)
f(x,t) = (e.*sin(w.*t)).*x.^2+(1-2.*e.*sin(w.*t)).*x;
dfdt = diff(f,t);
dfdx = diff(f,x);
dfdtdx = diff(diff(f,t),x)
dfdtdxdx = diff(diff(diff(f,t),x),x)
dfdxdx = diff(diff(f,x),x)
dfdxdxdx = diff(diff(diff(f,x),x),x)
%}
figure
subplot(121)
surface(s1,'edgecolor','none')
colorbar

subplot(122)
surface(lambda_0,'edgecolor','none')
colorbar


figure
subplot(121)
surface(l1,'edgecolor','none')
colorbar

subplot(122)
surface(lambda_1,'edgecolor','none')
colorbar


figure
subplot(121)
surface(a1,'edgecolor','none')
colorbar

subplot(122)
surface(lambda_2_first,'edgecolor','none')
colorbar

figure
subplot(121)
surface(a2,'edgecolor','none')
colorbar

subplot(122)
surface(lambda_2_second,'edgecolor','none')
colorbar



%{

%%% PLOT

A=0.1;
w=0.2*pi;
e=0.25;
t=0;
xdim=301;
ydim=151;
x=linspace(0,2,xdim);
y=linspace(0,1,ydim);
[x,y]=meshgrid(x,y);
f = (e.*sin(w.*t)).*x.^2+(1-2.*e.*sin(w.*t)).*x;
dfdt = (pi.*cos((pi.*t)./5).*x.^2)./20 - (pi.*cos((pi.*t)./5).*x)./10;
dfdx = (x.*sin((pi.*t)./5))./2 - sin((pi.*t)./5)./2 + 1;
dfdtdx = (pi.*x.*cos((pi.*t)./5))./10 - (pi.*cos((pi.*t)./5))./10
dfdtdxdx = (pi.*cos((pi.*t)./5))./10
dfdxdx = sin((pi.*t)./5)./2
dfdxdxdx = 0;


%%% S using arbitrary f(x,t), f(x,t) defined above
S_xx=-(A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdx)./2 - (A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdx)./2
S_xy= (A.*pi.^2.*sin(pi.*f).*sin(pi.*y))./2 - (A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.^2)./2 + (pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdx)./2
S_yy= (A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdx)./2 + (A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdx)./2

%%% B using arbitrary f(x,t), f(x,t) defined above
B_xx = (A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.^2 - pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdxdx).*(A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.^2 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdx) ...
 + (A.^2.*pi.^4.*cos(pi.*y).^2.*cos(pi.*f).^2.*dfdx.^2)./2 - (A.^2.*pi.^4.*cos(pi.*y).^2.*sin(pi.*f).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*sin(pi.*y).^2.*cos(pi.*f).^2.*dfdx.^2)./2 ...
 - (A.^2.*pi.^4.*sin(pi.*y).^2.*sin(pi.*f).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*cos(pi.*f).^2.*sin(pi.*y).^2.*dfdx.^2)./2 ...
 - (A.^2.*pi.^4.*sin(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2)./2 - (A.^2.*pi.^4.*sin(pi.*f).^2.*sin(pi.*y).^2.*dfdx.^2)./2 - (A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdtdx)./2 ...
 - (A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdtdx)./2 + (A.^2.*pi.^3.*cos(pi.*f).*sin(pi.*f).*cos(pi.*y).^2.*dfdxdx)./2 + (A.^2.*pi.^3.*cos(pi.*f).*sin(pi.*f).*sin(pi.*y).^2.*dfdxdx)./2 ...
 + (A.^2.*pi.^3.*cos(pi.*y).^2.*cos(pi.*f).*sin(pi.*f).*dfdxdx)./2 + (A.^2.*pi.^3.*sin(pi.*y).^2.*cos(pi.*f).*sin(pi.*f).*dfdxdx)./2 + (A.*pi.^3.*sin(pi.*f).*cos(pi.*y).*dfdx.*dfdt)./2 ...
 + (A.*pi.^3.*cos(pi.*y).*sin(pi.*f).*dfdx.*dfdt)./2 + A.^2.*pi.^4.*cos(pi.*f).*cos(pi.*y).^2.*cos(pi.*f).*dfdx.*dfdx;

B_xy = (pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdtdxdx)./2 + (A.*pi.*cos(pi.*y).*sin(pi.*f).*(A.*pi.^3.*sin(pi.*y).*cos(pi.*f).*dfdx.^3 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdxdx ...
 + 3.*A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.*dfdxdx))./2 + (A.*pi.^3.*cos(pi.*f).*sin(pi.*y).*dfdt)./2 - (A.*pi.^3.*sin(pi.*y).*cos(pi.*f).*dfdx.^2.*dfdt)./2 ...
 + (A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdx.*(A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.^2 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdx))./2 ...
 - A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdx.*(A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.^2 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdx) - A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.*dfdtdx ...
 - (A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdt.*dfdxdx)./2 - A.^2.*pi.^4.*sin(pi.*f).*cos(pi.*y).*sin(pi.*y).*cos(pi.*f).*dfdx - A.^2.*pi.^4.*cos(pi.*y).*sin(pi.*y).*cos(pi.*f).*sin(pi.*f).*dfdx.^3 ...
 + A.^2.*pi.^3.*cos(pi.*y).*sin(pi.*y).*cos(pi.*f).^2.*dfdx.*dfdxdx;
 
B_yx = (pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdtdxdx)./2 + (A.*pi.*sin(pi.*f).*cos(pi.*y).*(A.*pi.^3.*cos(pi.*f).*sin(pi.*y).*dfdx.^3 - pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdxdxdx ...
 + 3.*A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.*dfdxdx))./2 + (A.*pi.^3.*sin(pi.*y).*cos(pi.*f).*dfdt)./2 + (A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*(A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.^2 ...
 - pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdxdx).*dfdx)./2 - A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.*dfdtdx - (A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdt.*dfdxdx)./2 ...
 - A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdx.*(A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.^2 - pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdxdx) - (A.*pi.^3.*cos(pi.*f).*sin(pi.*y).*dfdx.^2.*dfdt)./2 ...
 - A.^2.*pi.^4.*cos(pi.*f).*sin(pi.*f).*cos(pi.*y).*sin(pi.*y).*dfdx.^3 - A.^2.*pi.^4.*cos(pi.*f).*cos(pi.*y).*sin(pi.*y).*sin(pi.*f).*dfdx + A.^2.*pi.^3.*cos(pi.*f).^2.*cos(pi.*y).*sin(pi.*y).*dfdx.*dfdxdx;

B_yy = (A.^2.*pi.^4.*cos(pi.*y).^2.*cos(pi.*f).^2.*dfdx.^2)./2 - (A.^2.*pi.^4.*sin(pi.*y).^2.*cos(pi.*f).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2)./2 ...
 + (A.*pi.*cos(pi.*y).*sin(pi.*f).*(A.*pi.^3.*cos(pi.*y).*sin(pi.*f).*dfdx.^2 - A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdxdx))./2 - (A.^2.*pi.^4.*cos(pi.*f).^2.*sin(pi.*y).^2.*dfdx.^2)./2 ...
 + (A.*pi.*sin(pi.*f).*cos(pi.*y).*(A.*pi.^3.*sin(pi.*f).*cos(pi.*y).*dfdx.^2 - A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdxdx))./2 - (A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*(A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.^2 ...
 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdx))./2 + (A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdtdx)./2 - (A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*(A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.^2 ...
 - pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdxdx))./2 + (A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdtdx)./2 + A.^2.*pi.^4.*sin(pi.*f).*sin(pi.*y).^2.*sin(pi.*f) - (A.*pi.^3.*sin(pi.*f).*cos(pi.*y).*dfdx.*dfdt)./2 ...
 - (A.*pi.^3.*cos(pi.*y).*sin(pi.*f).*dfdx.*dfdt)./2 + A.^2.*pi.^4.*cos(pi.*f).*cos(pi.*y).^2.*cos(pi.*f).*dfdx.*dfdx;


for i=1:ydim
    for j=1:xdim
        [V,D]=eig([S_xx(i,j),S_xy(i,j);S_xy(i,j),S_yy(i,j)]);
        if ~issorted(diag(D))
            [D,I] = sort(diag(D));
            V = V(:, I);
        end
        s1(i,j) = D(1,1);
        X0 = V(:,1);
        xi(i,j,:) = V(:,1);
        B=[B_xx(i,j),B_xy(i,j);B_xy(i,j),B_yy(i,j)];
        l1(i,j) = X0'*B*X0;
    end
end

figure
subplot(221)
surface(x,y,-s1,'edgecolor','none')
title('-s_1')
subplot(222)
surface(x,y,-s1.^2+0.5.*l1,'edgecolor','none')
title('correction')
subplot(223)
surface(x,y,xi(:,:,1),'edgecolor','none')
title('Xi_x')
subplot(224)
surface(x,y,xi(:,:,2),'edgecolor','none')
title('Xi_y')


figure
subplot(221)
surface(x,y,B_xx,'edgecolor','none')
title('B\_xx')
colorbar
subplot(222)
surface(x,y,B_xy,'edgecolor','none')
title('B\_xy')
colorbar
subplot(223)
surface(x,y,B_yx,'edgecolor','none')
title('B\_yx')
colorbar
subplot(224)
surface(x,y,B_yy,'edgecolor','none')
title('B\_yy')
colorbar
%
figure
quiver(x(1:15:end,1:15:end),y(1:15:end,1:15:end),xi(1:15:end,1:15:end,1),xi(1:15:end,1:15:end,2))%,'edgecolor','none')

figure
subplot(221)
surface(x,y,l1,'edgecolor','none')

%}