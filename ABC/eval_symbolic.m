close all
clear all
clc


%Compute Lambda terms for double gyre flow on domain [0 2]x[0 1]
dim = 96
xx = linspace(0,2*pi,dim);
[xx,yy,zz]=meshgrid(xx,xx,xx);

tt = 0;
for i = 1:dim
    i
    for j = 1:dim
        for k = 1:dim
            x=xx(i,j,k);
            y=yy(i,j,k);
            z=zz(i,j,k);
            
            SS = [0, (2^(1/2)*cos(x))/2 - sin(y)/2, (3^(1/2)*cos(z))/2 - (2^(1/2)*sin(x))/2;
            (2^(1/2)*cos(x))/2 - sin(y)/2, 0, cos(y)/2 - (3^(1/2)*sin(z))/2;
            (3^(1/2)*cos(z))/2 - (2^(1/2)*sin(x))/2, cos(y)/2 - (3^(1/2)*sin(z))/2, 0];            

            BB = [2*cos(x)^2 + 2*sin(x)^2 - 2^(1/2)*cos(x)*sin(y) - 2^(1/2)*3^(1/2)*cos(z)*sin(x), (3^(1/2)*cos(y)*cos(z))/2 - (cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/2 - 2^(1/2)*cos(y)*sin(x) - (2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)))/2 + (2^(1/2)*3^(1/2)*sin(x)*sin(z))/2, (2^(1/2)*cos(x)*cos(y))/2 - (2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z)))/2 - (3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)))/2 + (3^(1/2)*sin(y)*sin(z))/2 - 2^(1/2)*3^(1/2)*cos(x)*sin(z);
            (3^(1/2)*cos(y)*cos(z))/2 - (cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/2 - 2^(1/2)*cos(y)*sin(x) - (2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)))/2 + (2^(1/2)*3^(1/2)*sin(x)*sin(z))/2, cos(y)^2 + sin(y)^2 - 2^(1/2)*cos(x)*sin(y) - 3^(1/2)*cos(y)*sin(z), (2^(1/2)*sin(x)*sin(y))/2 - (3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)))/2 - 3^(1/2)*cos(z)*sin(y) - (sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/2 + (2^(1/2)*3^(1/2)*cos(x)*cos(z))/2;
            (2^(1/2)*cos(x)*cos(y))/2 - (2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z)))/2 - (3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)))/2 + (3^(1/2)*sin(y)*sin(z))/2 - 2^(1/2)*3^(1/2)*cos(x)*sin(z), (2^(1/2)*sin(x)*sin(y))/2 - (3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)))/2 - 3^(1/2)*cos(z)*sin(y) - (sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/2 + (2^(1/2)*3^(1/2)*cos(x)*cos(z))/2, 3*cos(z)^2 + 3*sin(z)^2 - 3^(1/2)*cos(y)*sin(z) - 2^(1/2)*3^(1/2)*cos(z)*sin(x)]; 

            QQ = [(2*(2^(1/2)*sin(x)*sin(y) - 2^(1/2)*3^(1/2)*cos(x)*cos(z))*(cos(y) + 3^(1/2)*sin(z)))/3 - 2*2^(1/2)*sin(x)*(2^(1/2)*cos(x)*cos(y) - 2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z))) + (2*2^(1/2)*sin(x)*(3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*sin(y)*sin(z)))/3 - (2*2^(1/2)*cos(x)*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(y)*cos(z)))/3 - 2*2^(1/2)*cos(x)*(2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*sin(x)*sin(z)) - (2*2^(1/2)*cos(x)*cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 + (2*2^(1/2)*3^(1/2)*sin(x)*sin(z)*(sin(y) + 2^(1/2)*cos(x)))/3, (4*sin(y)*(2^(1/2)*cos(x)*sin(y) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 - ((2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*cos(x)*sin(z))*(cos(y) + 3^(1/2)*sin(z)))/3 + cos(y)*(2^(1/2)*cos(x)*cos(y) - 2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z))) - (cos(y)*(3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*sin(y)*sin(z)))/3 + ((3^(1/2)*cos(z) + 2^(1/2)*sin(x))*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(z)*sin(y)))/3 - (4*2^(1/2)*cos(x)*(2^(1/2)*cos(x)*sin(y) + 3^(1/2)*cos(y)*sin(z)))/3 + (2^(1/2)*sin(x)*(3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*cos(z)))/3 + 2^(1/2)*sin(x)*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*sin(x)*sin(y)) + (2^(1/2)*sin(x)*sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 - (2^(1/2)*cos(x)*cos(y)*(cos(y) + 3^(1/2)*sin(z)))/3, (3^(1/2)*sin(z)*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(y)*cos(z)))/3 - ((2^(1/2)*cos(y)*sin(x) - 2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)))*(cos(y) + 3^(1/2)*sin(z)))/3 - 2^(1/2)*cos(x)*(3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*cos(z)) - ((3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*cos(z)*sin(y))*(sin(y) + 2^(1/2)*cos(x)))/3 - (4*3^(1/2)*cos(z)*(2^(1/2)*cos(x)*sin(y) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 + (4*2^(1/2)*sin(x)*(3^(1/2)*cos(y)*sin(z) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 - (2^(1/2)*cos(x)*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*sin(x)*sin(y)))/3 + 3^(1/2)*sin(z)*(2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*sin(x)*sin(z)) - (2^(1/2)*3^(1/2)*cos(x)*cos(z)*(sin(y) + 2^(1/2)*cos(x)))/3 + (2^(1/2)*3^(1/2)*sin(x)*sin(z)*(cos(y) + 3^(1/2)*sin(z)))/3;
            (4*sin(y)*(2^(1/2)*cos(x)*sin(y) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 - ((2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*cos(x)*sin(z))*(cos(y) + 3^(1/2)*sin(z)))/3 + cos(y)*(2^(1/2)*cos(x)*cos(y) - 2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z))) - (cos(y)*(3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*sin(y)*sin(z)))/3 + ((3^(1/2)*cos(z) + 2^(1/2)*sin(x))*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(z)*sin(y)))/3 - (4*2^(1/2)*cos(x)*(2^(1/2)*cos(x)*sin(y) + 3^(1/2)*cos(y)*sin(z)))/3 + (2^(1/2)*sin(x)*(3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*cos(z)))/3 + 2^(1/2)*sin(x)*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*sin(x)*sin(y)) + (2^(1/2)*sin(x)*sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 - (2^(1/2)*cos(x)*cos(y)*(cos(y) + 3^(1/2)*sin(z)))/3, 2*sin(y)*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(y)*cos(z)) - 2*cos(y)*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*sin(x)*sin(y)) + (2*sin(y)*(2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*sin(x)*sin(z)))/3 - (2*(2^(1/2)*cos(x)*cos(y) - 3^(1/2)*sin(y)*sin(z))*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 - (2*cos(y)*(3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*cos(z)))/3 - (2*3^(1/2)*cos(y)*cos(z)*(sin(y) + 2^(1/2)*cos(x)))/3 + (2*2^(1/2)*sin(x)*sin(y)*(cos(y) + 3^(1/2)*sin(z)))/3,((3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*sin(z))*(sin(y) + 2^(1/2)*cos(x)))/3 - (4*cos(y)*(3^(1/2)*cos(y)*sin(z) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 - (sin(y)*(2^(1/2)*cos(x)*cos(y) - 2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z))))/3 - ((3^(1/2)*cos(z) + 2^(1/2)*sin(x))*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*cos(y)*sin(x)))/3 + sin(y)*(3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*sin(y)*sin(z)) + (4*3^(1/2)*sin(z)*(2^(1/2)*cos(x)*sin(y) + 3^(1/2)*cos(y)*sin(z)))/3 - 3^(1/2)*cos(z)*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(y)*cos(z)) - (3^(1/2)*cos(z)*(2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*sin(x)*sin(z)))/3 - (3^(1/2)*cos(y)*cos(z)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 + (3^(1/2)*sin(y)*sin(z)*(sin(y) + 2^(1/2)*cos(x)))/3;
            (3^(1/2)*sin(z)*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(y)*cos(z)))/3 - ((2^(1/2)*cos(y)*sin(x) - 2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)))*(cos(y) + 3^(1/2)*sin(z)))/3 - 2^(1/2)*cos(x)*(3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*cos(z)) - ((3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*cos(z)*sin(y))*(sin(y) + 2^(1/2)*cos(x)))/3 - (4*3^(1/2)*cos(z)*(2^(1/2)*cos(x)*sin(y) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 + (4*2^(1/2)*sin(x)*(3^(1/2)*cos(y)*sin(z) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 - (2^(1/2)*cos(x)*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*sin(x)*sin(y)))/3 + 3^(1/2)*sin(z)*(2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*sin(x)*sin(z)) - (2^(1/2)*3^(1/2)*cos(x)*cos(z)*(sin(y) + 2^(1/2)*cos(x)))/3 + (2^(1/2)*3^(1/2)*sin(x)*sin(z)*(cos(y) + 3^(1/2)*sin(z)))/3, ((3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*sin(z))*(sin(y) + 2^(1/2)*cos(x)))/3 - (4*cos(y)*(3^(1/2)*cos(y)*sin(z) + 2^(1/2)*3^(1/2)*cos(z)*sin(x)))/3 - (sin(y)*(2^(1/2)*cos(x)*cos(y) - 2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z))))/3 - ((3^(1/2)*cos(z) + 2^(1/2)*sin(x))*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*cos(y)*sin(x)))/3 + sin(y)*(3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*sin(y)*sin(z)) + (4*3^(1/2)*sin(z)*(2^(1/2)*cos(x)*sin(y) + 3^(1/2)*cos(y)*sin(z)))/3 - 3^(1/2)*cos(z)*(cos(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 3^(1/2)*cos(y)*cos(z)) - (3^(1/2)*cos(z)*(2^(1/2)*sin(x)*(cos(y) + 3^(1/2)*sin(z)) - 2^(1/2)*3^(1/2)*sin(x)*sin(z)))/3 - (3^(1/2)*cos(y)*cos(z)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 + (3^(1/2)*sin(y)*sin(z)*(sin(y) + 2^(1/2)*cos(x)))/3, 2*3^(1/2)*sin(z)*(3^(1/2)*cos(z)*(sin(y) + 2^(1/2)*cos(x)) - 2^(1/2)*3^(1/2)*cos(x)*cos(z)) - 2*3^(1/2)*cos(z)*(3^(1/2)*sin(z)*(sin(y) + 2^(1/2)*cos(x)) - 3^(1/2)*sin(y)*sin(z)) - (2*(3^(1/2)*cos(y)*cos(z) - 2^(1/2)*3^(1/2)*sin(x)*sin(z))*(sin(y) + 2^(1/2)*cos(x)))/3 + (2*3^(1/2)*sin(z)*(sin(y)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)) - 2^(1/2)*sin(x)*sin(y)))/3 + (2*3^(1/2)*cos(z)*(2^(1/2)*cos(x)*cos(y) - 2^(1/2)*cos(x)*(cos(y) + 3^(1/2)*sin(z))))/3 + (2*3^(1/2)*sin(y)*sin(z)*(3^(1/2)*cos(z) + 2^(1/2)*sin(x)))/3 - (2*2^(1/2)*3^(1/2)*cos(x)*cos(z)*(cos(y) + 3^(1/2)*sin(z)))/3];
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