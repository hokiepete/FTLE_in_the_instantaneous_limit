close all
clear all
clc

syms u(x,y,t) v(x,y,t) f(x,t)
assume(x,'real')
assume(y,'real')
assume(t,'real')

A=0.1;
w=0.2*pi;
e=0.25;

f(x,t) = (e.*sin(w.*t)).*x.^2+(1-2.*e.*sin(w.*t)).*x;
%Velocity field
u(x,y,t) =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
v(x,y,t) = pi.*A.*cos(pi.*f).*sin(y.*pi).*diff(f,x);

%First Derivatives of velocity
dudt = diff(u,t);
dudx = diff(u,x);
dudy = diff(u,y);

dvdt = diff(v,t);
dvdx = diff(v,x);
dvdy = diff(v,y);

%Total acceleration field
au = dudt + u*dudx + v*dudy;
av = dvdt + u*dvdx + v*dvdy;

%Derivatives of acceleration
daudt = diff(au,t);
daudx = diff(au,x);
daudy = diff(au,y);

davdt = diff(av,t);
davdx = diff(av,x);
davdy = diff(av,y);

%Total jerk field
ju = daudt+u.*daudx+v.*daudy;
jv = davdt+u.*davdx+v.*davdy;

%Derivatives of jerk
djudt = diff(ju,t);
djudx = diff(ju,x);
djudy = diff(ju,y);

djvdt = diff(jv,t);
djvdx = diff(jv,x);
djvdy = diff(jv,y);

%Gradient matrices 
grad_v = [dudx,dudy;dvdx,dvdy];
grad_a = [daudx,daudy;davdx,davdy];
grad_j = [djudx,djudy;djvdx,djvdy];

%Rivlin-Ericksen tensors
S = 0.5*(grad_v+grad_v');
B = 0.5*(grad_a+grad_a')+grad_v'*grad_v;
Q = 0.5*(grad_j + grad_j')+(grad_v'*grad_a+grad_a'*grad_v);

%Rotation matrix
R = [0,-1;1,0];

%Compute Lambda terms for double gyre flow on domain [0 2]x[0 1]
xdim = 301
ydim =151
xx = linspace(0,2,xdim);
yy = linspace(0,1,ydim);
[xx,yy]=meshgrid(xx,yy);
tt = 0;
for i = 1:ydim
    for j = 1:xdim
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