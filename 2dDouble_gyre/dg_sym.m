close all
clear all
clc

%

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

daudx = diff(au,x);
daudy = diff(au,y);
davdx = diff(av,x);
davdy = diff(av,y);

grad_v = [dudx,dudy;dvdx,dvdy];
grad_a = [daudx,daudy;davdx,davdy];

S = 0.5*(grad_v+grad_v');
B = 0.5*(grad_a+grad_a')+grad_v'*grad_v;

%}

%{

%%% Known f(x,t)

syms u(x,y,t) v(x,y,t) f(x,t) e w A
assume(A,'real')
assume(e,'real')
assume(w,'real')
assume(x,'real')
assume(y,'real')
assume(t,'real')


f(x,t) = (e.*sin(w.*t)).*x.^2+(1-2.*e.*sin(w.*t)).*x;
dfdx = diff(f,x);
u(x,y,t) =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
v(x,y,t) = pi.*A.*cos(pi.*f).*sin(y.*pi).*dfdx;
dudt = diff(u,t);
dudx = diff(u,x);
dudy = diff(u,y);

dvdt = diff(v,t);
dvdx = diff(v,x);
dvdy = diff(v,y);

au = dudt + u*dudx + v*dudy;
av = dvdt + u*dvdx + v*dvdy;

daudx = diff(au,x);
daudy = diff(au,y);
davdx = diff(av,x);
davdy = diff(av,y);

grad_v = [dudx,dudy;dvdx,dvdy];
grad_a = [daudx,daudy;davdx,davdy];

S = 0.5*(grad_v+grad_v');
B = 0.5*(grad_a+grad_a')+grad_v'*grad_v;
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


%{

%%% PLOT

A=0.1;
w=0.2*pi;
e=0.25;
t=0;
x=linspace(0,2,301);
y=linspace(0,1,151);
[x,y]=meshgrid(x,y);
f = (e.*sin(w.*t)).*x.^2+(1-2.*e.*sin(w.*t)).*x;
dfdt = (pi.*cos((pi.*t)./5).*x.^2)./20 - (pi.*cos((pi.*t)./5).*x)./10;
dfdx = (x.*sin((pi.*t)./5))./2 - sin((pi.*t)./5)./2 + 1;
dfdtdx = (pi.*x.*cos((pi.*t)./5))./10 - (pi.*cos((pi.*t)./5))./10
dfdtdxdx = (pi.*cos((pi.*t)./5))./10
dfdxdx = sin((pi.*t)./5)./2
dfdxdxdx = 0;

%%% B using arbitrary f(x,t), f(x,t) defined above
B_xx = (A.*pi.^2.*sin(pi.*f).*sin(pi.*y).*dfdx.^2 - pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdxdx).*(A.*pi.^2.*sin(pi.*y).*sin(pi.*f).*dfdx.^2 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdx) ...
 + (A.^2.*pi.^4.*cos(pi.*y).^2.*cos(pi.*f).^2.*dfdx.^2)./2 - (A.^2.*pi.^4.*cos(pi.*y).^2.*sin(pi.*f).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*sin(pi.*y).^2.*cos(pi.*f).^2.*dfdx.^2)./2 ...
 - (A.^2.*pi.^4.*sin(pi.*y).^2.*sin(pi.*f).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2)./2 + (A.^2.*pi.^4.*cos(pi.*f).^2.*sin(pi.*y).^2.*dfdx.^2)./2 ...
 - (A.^2.*pi.^4.*sin(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2)./2 - (A.^2.*pi.^4.*sin(pi.*f).^2.*sin(pi.*y).^2.*dfdx.^2)./2 - (A.*pi.^2.*cos(pi.*y).*cos(pi.*f).*dfdtdx)./2 ...
 - (A.*pi.^2.*cos(pi.*f).*cos(pi.*y).*dfdtdx)./2 + (A.^2.*pi.^3.*cos(pi.*f).*sin(pi.*f).*cos(pi.*y).^2.*dfdxdx)./2 + (A.^2.*pi.^3.*cos(pi.*f).*sin(pi.*f).*sin(pi.*y).^2.*dfdxdx)./2 ...
 + (A.^2.*pi.^3.*cos(pi.*y).^2.*cos(pi.*f).*sin(pi.*f).*dfdxdx)./2 + (A.^2.*pi.^3.*sin(pi.*y).^2.*cos(pi.*f).*sin(pi.*f).*dfdxdx)./2 + (A.*pi.^3.*sin(pi.*f).*cos(pi.*y).*dfdx.*dfdt)./2 ...
 + (A.*pi.^3.*cos(pi.*y).*sin(pi.*f).*dfdx.*dfdt)./2 + A.^2.*pi.^4.*cos(pi.*f).*cos(pi.*y).^2.*cos(pi.*f).*dfdx.*dfdx;

B_xy = 3(pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdtdxdx)./2 + (A.*pi.*cos(pi.*y).*sin(pi.*f).*(A.*pi.^3.*sin(pi.*y).*cos(pi.*f).*dfdx.^3 - pi.*A.*sin(pi.*y).*cos(pi.*f).*dfdxdxdx ...
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
 
subplot(221)
surface(x,y,B_xx,'edgecolor','none')
title('B\_xx')
subplot(222)
surface(x,y,B_xy,'edgecolor','none')
title('B\_xy')
subplot(223)
surface(x,y,B_yx,'edgecolor','none')
title('B\_yx')
subplot(224)
surface(x,y,B_yy,'edgecolor','none')
title('B\_yy')
%}

