close all
clear all
clc

t=0.5*pi
x=linspace(0,2,301);
y=linspace(0,1,151);

[x,y]=meshgrid(x,y);
A = 0.1;
w = 0.2.*pi;
e = 0.25;
a = e.*sin(w.*t);
b = 1-2.*e.*sin(w.*t);
f = a.*x.^2+b.*x;
dfdx = 2.*a.*x+b;
dfdy = zeros(size(x));
dfdt = e.*w.*cos(w.*t).*x.^2-2.*e.*w.*cos(w.*t).*x;
dfdtdx= 2.*e.*w.*cos(w.*t).*x-2.*e.*w.*cos(w.*t);
dfdxdx = 2.*a+zeros(size(x));
dfdtdxdx = 2.*e.*w.*cos(w.*t)+zeros(size(x));

u =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
v = pi.*A.*cos(pi.*f).*sin(y.*pi).*dfdx;

dudt = -pi.^2.*A.*cos(pi.*f).*cos(pi.*y).*dfdt+0.5.*pi.^3.*A.^2.*sin(2.*pi.*f).*dfdx;
 
dvdt = pi.^2.*A.*(-sin(pi.*f).*sin(pi.*y).*dfdx.*dfdt+1./pi.*cos(pi.*f).*sin(pi.*y).*dfdtdx)...
    +0.5.*pi.^3.*A.^2.*sin(2.*pi.*y).*(sin(pi.*f).^2.*dfdx+cos(pi.*f).^2.*dfdx.^2-1./(2.*pi).*sin(2.*pi.*f).*dfdxdx);
        
B_xx = pi.^2.*A.*cos(pi.*y).*(pi.*sin(pi.*f).*dfdx.*dfdt-cos(pi.*f).*dfdtdx)...
+pi.^4.*A.^2.*(dfdx.^2.*(cos(2.*pi.*f)+cos(pi.*f).^2.*cos(pi.*y).^2+sin(pi.*f).^2.*sin(pi.*y).^2)...
+1./pi.*sin(2.*pi.*f).*dfdxdx.*(1./2-dfdx)+1./(pi.^2).*cos(pi.*f).^2.*sin(pi.*y).^2.*dfdxdx.^2);

B_xy = 0.5.*pi.^2.*A.*sin(pi.*y).*(sin(pi.*f).*(pi.*dfdt-dfdxdx.*dfdt-2.*dfdx.*dfdtdx)...
+cos(pi.*f).*(-pi.*dfdx.^2.*dfdt+1./pi.*dfdtdxdx))...
+0.25.*pi.^4.*A.^2.*sin(2.*pi.*y).*(-sin(2.*pi.*f).*dfdx.*(1+dfdx.^2)...
+1./pi.*dfdxdx.*(sin(pi.*f).^2+4.*cos(pi.*f).^2.*dfdx-cos(2.*pi.*f).*dfdx));

B_yy = pi.^2.*A.*cos(pi.*y).*(-pi.*sin(pi.*f).*dfdx.*dfdt+cos(pi.*f).*dfdtdx)...
+pi.^4.*A.^2.*(cos(2.*pi.*f).*(sin(pi.*f).^2.*dfdx+cos(pi.*f).^2.*dfdx.^2-1./(2.*pi).*sin(2.*pi.*f).*dfdxdx)...
+sin(pi.*f).^2.*sin(pi.*y).^2+cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2);

alpha = sin(pi.*f).*sin(pi.*y).*(1-dfdx)+1./pi.*cos(pi.*f).*sin(pi.*y).*dfdxdx;
beta = cos(pi.*f).*cos(pi.*y).*dfdx;

s1_bar = -0.5.*sqrt(alpha.^2+4.*beta.^2);
N = sqrt(0.25.*alpha.^2+(s1_bar+beta).^2);

xi_x=1./N.*0.5.*alpha;
xi_y=1./N.*s1_bar+beta;

s1=-0.5.*pi.^2.*A.*sqrt(cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2 ...
    +0.25.*sin(pi.*f).^2.*sin(pi.*y).^2.*(1-dfdx).^2);

correction = s1.^2-0.5.*(B_xx.*xi_x.^2+2.*B_xy.*xi_x.*xi_y+B_yy.*xi_y.^2);

subplot(621)
surface(x,y,s1,'edgecolor','none')
title('s1')
subplot(622)
surface(x,y,correction,'edgecolor','none')
title('correction')
subplot(623)
surface(x,y,xi_x,'edgecolor','none')
title('xi\_x')
subplot(624)
surface(x,y,xi_y,'edgecolor','none')
title('xi\_y')
subplot(625)
surface(x,y,B_xx,'edgecolor','none')
title('B\_xx')
subplot(626)
surface(x,y,B_xy,'edgecolor','none')
title('B\_xy')
subplot(627)
surface(x,y,B_yy,'edgecolor','none')
title('B\_yy')
subplot(628)
surface(x,y,N,'edgecolor','none')
title('N')
subplot(629)
surface(x,y,alpha,'edgecolor','none')
title('alpha')
subplot(6,2,10)
surface(x,y,beta,'edgecolor','none')
title('beta')
subplot(6,2,11)
surface(x,y,s1_bar,'edgecolor','none')
title('s1\_bar')





