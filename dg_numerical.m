close all
clear all
clc
t=0
x=linspace(0,2,301)
dx=x(2)-x(1);
y=linspace(0,1,151)
dy=y(2)-y(1);
[x,y]=meshgrid(x,y)
A = 0.1;
w = 0.2.*pi;
e = 0.25;
a = e.*sin(w.*t);
b = 1-2.*e.*sin(w.*t);
f = a.*x.^2+b.*x;
dfdx = 2.*a.*x+b;    
u =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
v = pi.*A.*cos(pi.*f).*sin(y.*pi).*dfdx;

[dudx,dudy] = gradient(u,dx,dy);
[dvdx,dvdy] = gradient(v,dx,dy);

for i = 1:151
    i
    for j=1:301
        %utemp =[u(i,j);v(i,j);w(i,j)];
        gradU = [dudx(i,j),dudy(i,j);dvdx(i,j),dvdy(i,j)];
        S = 0.5*(gradU + gradU');
        %s1(i,j) = min(eig(S));

        %s3(i,j) = max(eig(S));
        [V,D] = eig(S);
        if ~issorted(diag(D))
            [D,I] = sort(diag(D));
            V = V(:, I);
        end
        s1(i,j) = D(1,1);
        X1(i,j,:) = V(:,1);
        sn(i,j) = D(end,end);
        Xn(i,j,:) = V(:,end);
        %s1(i,j) = median(eig(S));
    end
end

surface(x,y,s1,'edgecolor','none')
