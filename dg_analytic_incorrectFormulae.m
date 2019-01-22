t=0
x=linspace(0,2,301)
y=linspace(0,1,151)

[x,y]=meshgrid(x,y)
A = 0.1
w = 0.2.*pi
e = 0.25
a = e.*sin(w.*t)
b = 1-2.*e.*sin(w.*t)
f = a.*x.^2+b.*x
dfdx = 2.*a.*x+b    
dfdy = 0
dfdt = e.*w.*cos(w.*t).*x.^2-2.*e.*w.*cos(w.*t).*x
dfdtdx= 2.*e.*w.*cos(w.*t).*x-2.*e.*w.*cos(w.*t)
dfdxdx = 2.*a
dfdtdxdx = 2.*e.*w.*cos(w.*t)

u =-pi.*A.*sin(pi.*f).*cos(y.*pi)    
v = pi.*A.*cos(pi.*f).*sin(y.*pi).*dfdx

dudt = -pi.^2.*A.*cos(pi.*f).*cos(pi.*y).*dfdt+0.5.*pi.^3.*A.^2.*sin(2.*pi.*f).*dfdx
 
dvdt = -pi.^2.*A.*sin(pi.*f).*sin(pi.*y).*dfdx.*dfdt+pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdtdx...
    +0.5.*pi.^3.*A.^2.*sin(2.*pi.*y).*dfdx.*(sin(pi.*f).^2+cos(pi.*f).^2.*dfdx)


B_xx = pi.^3.*A.*sin(pi.*f).*cos(pi.*y).*dfdx.*dfdt-pi.^2.*A.*cos(pi.*f).*cos(pi.*y).*dfdtdx+0.5.*pi.^3.*A.^2.*sin(2.*pi.*y).*dfdxdx...
    +pi.^4.*A.^2.*dfdx.^2.*(cos(pi.*f).^2.*cos(pi.*y).^2+sin(pi.*f).^2.*sin(pi.*y).^2+sin(2.*pi.*f))

    
B_xy = 0.5.*pi.^3.*A.*cos(pi.*f).*sin(pi.*y).*dfdt.*(1-dfdx.^2)...
    -0.5.*pi.^2.*A.*sin(pi.*f).*sin(pi.*y).*dfdtdx.*(dfdt+2.*dfdx)...
    +0.5.*pi.*A.*cos(pi.*f).*sin(pi.*y).*dfdtdxdx...
    +0.25.*pi.^3.*A.^2.*sin(2.*pi.*y).*dfdxdx.*(sin(pi.*f).^2+cos(pi.*f).^2.*dfdx)...
    +0.25.*pi.^3.*A.^2.*sin(2.*pi.*y).*dfdx.*...
    (-pi.*sin(2.*pi.*f).*(1+dfdx.^2)+cos(pi.*f).^2.*dfdxdx)...
    +0.5.*pi.^4.*A.^2.*cos(2.*pi.*y).*dfdx.*...
    (sin(pi.*f).^2+cos(pi.*f).^2.*dfdx)
     
B_yy = -pi.^3.*A.*sin(pi.*f).*cos(pi.*y).*dfdx.*dfdt...
    +pi.^2.*A.*cos(pi.*f).*cos(pi.*y).*dfdtdx...
    +pi.^4.*A.^2.*cos(2.*pi.*y).*dfdx.*(sin(pi.*f).^2+cos(pi.*f).^2.*dfdx)...
    +pi.^4.*A.^2.*(cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2 ...
    +sin(pi.*f).^2.*sin(pi.*y).^2)

alpha = sin(pi.*f).*sin(pi.*y).*(1-dfdx)
beta = cos(pi.*f).*cos(pi.*y).*dfdx

s1_bar = -0.5.*sqrt(0.25.*alpha.^2+beta.^2)
N = sqrt(0.25.*alpha.^2+(s1_bar+beta).^2)

xi_x=1./N.*0.5.*alpha
xi_y=1./N.*s1_bar+beta

s1=-0.5.*pi.^2.*A.*sqrt(cos(pi.*f).^2.*cos(pi.*y).^2.*dfdx.^2 ...
    +0.25.*sin(pi.*f).^2.*sin(pi.*y).^2.*(1-dfdx).^2)


correction = s1.^2-0.5.*(B_xx.*xi_x.^2+2.*B_xy.*xi_x.*xi_y+B_yy.*xi_y.^2)

