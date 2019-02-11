

%%% Start with velocity field u and v

%velocity derivatives
[dudx,dudy,dudt] = gradient(u,dx,dy,dt);
[dvdx,dvdy,dvdt] = gradient(v,dx,dy,dt);

au = dudt+u.*dudx+v.*dudy;
av = dvdt+u.*dvdx+v.*dvdy;

%acceleration derivatives
[daudx,daudy,daudt] = gradient(au,dx,dy,dt);
[davdx,davdy,davdt] = gradient(av,dx,dy,dt);

ju = daudt+u.*daudx+v.*daudy;
jv = davdt+u.*davdx+v.*davdy;

%jerk derivatives
[djudx,djudy,djudt] = gradient(ju,dx,dy,dt);
[djvdx,djvdy,djvdt] = gradient(jv,dx,dy,dt);

R = [0,-1;1,0];

%Loop through each grid point
for t =1:tdim
    for i = 1:ydim
        for j = 1:xdim
            %check for NaNs
            if ~isnan(ju(i,j,t))&&~isnan(jv(i,j,t))
                %gradient of velocity
                Grad_v = [dudx(i,j,t), dudy(i,j,t);dvdx(i,j,t), dvdy(i,j,t)];
                
                %gradient of accelerations 
                Grad_a = [daudx(i,j,t), daudy(i,j,t); davdx(i,j,t), davdy(i,j,t)];
                
                %gradient of jerk
                Grad_j = [djudx(i,j,t), djudy(i,j,t); djvdx(i,j,t), djvdy(i,j,t)];
                
                %Tensors
                S = 0.5*(Grad_v + Grad_v');
                B = 0.5*(Grad_a + Grad_a')+(Grad_v'*Grad_v);
                Q = 0.3*(Grad_j + Grad_j')+(Grad_v'*Grad_a+Grad_a'*Grad_v);
                
                %Calculate s1, lambda1, and lambda2 
                [V,D] = eig(S);
                if ~issorted(diag(D))
                    [D,I] = sort(diag(D));
                    V = V(:, I);
                end
                s1(i,j,t) = D(1,1);
                X0 = V(:,1);
                lambda_1(i,j,t) = X0'*B*X0;
                
                %First lambda_2 method calculating Xi_1
                X1 = pinv((S - s1(i,j,t)*eye(size(B))))*(-((B - l1(i,j,t)*eye(size(B)))*X0));
                
                lambda_2_first(i,j,t) = X0'*Q*X0 + X0'*B*X1 - X0'*S*X1;
                
                %Second lambda_2 method bypassing Xi_1
                mu = X0'*R'*(S-s1(i,j,t)*eye(size(S)))*R*X0;
                d = X0'*R'*B*X0;
                lambda_2_second(i,j,t) = X0'*Q*X0 - d.^2/m;
               
                
            else
                s1(i,j,t) =nan;
                lambda_1(i,j,t) = nan;
                lambda_2_first(i,j,t) = nan;
                lambda_2_second(i,j,t) = nan;
            end
        end
    end
end
