close all
clear all
clc

syms u(x,y,z,t) v(x,y,z,t) w(x,y,z,t)
assume(x,'real')
assume(y,'real')
assume(z,'real')
assume(t,'real')

ABC_Amplitude=0.0;
Ap = sqrt(3);
Bp = sqrt(2);
%Velocity field
u(x,y,z,t) = (Ap+ABC_Amplitude*sin(pi*t))*sin(z) + cos(y);    
v(x,y,z,t) = Bp*sin(x) + (Ap+ABC_Amplitude*sin(pi*t))*cos(z);
w(x,y,z,t) = sin(y) + Bp*cos(x);

%First Derivatives of velocity
dudt = diff(u,t);
dudx = diff(u,x);
dudy = diff(u,y);
dudz = diff(u,z);

dvdt = diff(v,t);
dvdx = diff(v,x);
dvdy = diff(v,y);
dvdz = diff(v,z);

dwdt = diff(w,t);
dwdx = diff(w,x);
dwdy = diff(w,y);
dwdz = diff(w,z);

%Total acceleration field
au = dudt + u*dudx + v*dudy + w*dudz;
av = dvdt + u*dvdx + v*dvdy + w*dvdz;
aw = dwdt + u*dwdx + v*dwdy + w*dwdz;

%Derivatives of acceleration
daudt = diff(au,t);
daudx = diff(au,x);
daudy = diff(au,y);
daudz = diff(au,z);

davdt = diff(av,t);
davdx = diff(av,x);
davdy = diff(av,y);
davdz = diff(av,z);

dawdt = diff(aw,t);
dawdx = diff(aw,x);
dawdy = diff(aw,y);
dawdz = diff(aw,z);

%Total jerk field
ju = daudt + u*daudx + v*daudy + w*daudz;
jv = davdt + u*davdx + v*davdy + w*davdz;
jw = dawdt + u*dawdx + v*dawdy + w*dawdz;

%Derivatives of jerk
djudt = diff(ju,t);
djudx = diff(ju,x);
djudy = diff(ju,y);
djudz = diff(ju,z);

djvdt = diff(jv,t);
djvdx = diff(jv,x);
djvdy = diff(jv,y);
djvdz = diff(jv,z);

djwdt = diff(jw,t);
djwdx = diff(jw,x);
djwdy = diff(jw,y);
djwdz = diff(jw,z);

%Gradient matrices 
grad_v = [dudx,dudy,dudz;dvdx,dvdy,dvdz;dwdx,dwdy,dwdz];
grad_a = [daudx,daudy,daudz;davdx,davdy,davdz;dawdx,dawdy,dawdz];
grad_j = [djudx,djudy,djudz;djvdx,djvdy,djvdz;djwdx,djwdy,djwdz];

%Rivlin-Ericksen tensors
S = 0.5*(grad_v+grad_v');
B = 0.5*(grad_a+grad_a')+grad_v'*grad_v;
Q = 1./3.*(grad_j + grad_j')+(grad_v'*grad_a+grad_a'*grad_v);


%{
%Compute Lambda terms for double gyre flow on domain [0 2]x[0 1]
dim = 192
x = linspace(0,2*pi,dim);
[x,y,z]=meshgrid(x,x,x);

tt = 0;
for i = 1:dim
    i
    for j = 1:dim
        for k = 1:dim
            x=xx(i,j,k);
            x=xx(i,j,k);
            x=xx(i,j,k);
            SS = reshape(double(subs(S,[x,y,z,t],[xx(i,j,k),yy(i,j,k),zz(i,j,k),tt])),3,3);
            BB = reshape(double(subs(B,[x,y,z,t],[xx(i,j,k),yy(i,j,k),zz(i,j,k),tt])),3,3);
            QQ = reshape(double(subs(Q,[x,y,z,t],[xx(i,j,k),yy(i,j,k),zz(i,j,k),tt])),3,3);
            [V D] = eig(SS);
            if ~issorted(diag(D))
                [D,I] = sort(diag(D));
                V = V(:, I);
            end
            lambda_0(i,j,k) = D(1,1);
            X0 = V(:,1);
            lambda_1(i,j,k) = X0'*BB*X0;

            %First lambda_2 method calculating Xi_1
            %X1 = -((B-l1(i,j,t)*eye(size(B)))*X0)\(S-s1(i,j,t)*eye(size(B)));
            X1 = pinv((SS-lambda_0(i,j,k)*eye(size(BB))))*(-(BB-lambda_1(i,j,k)*eye(size(BB)))*X0);
            %X1=X1';
            %{
            if sum(X1)~=0
                X1=X1'/norm(X1);
            else
                X1=X1';
            end
            %}        
            %lambda_2_first(i,j) = X0'*QQ*X0 + X0'*BB*X1 - X0'*SS*X1;
            lambda_2(i,j,k) = X0'*QQ*X0 + X0'*BB*X1 - lambda_1(i,j,k).*X0'*X1;
        end
    end
end
%{
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
%}

save analytic_lambda_terms lambda_0 lambda_1 lambda_2
%}