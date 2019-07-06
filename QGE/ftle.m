close all
clear all
clc
load flow_map_qge_v2
tt=abs(time-time(1));
dx = x(1,2)-x(1,1)
dy = y(2,1)-y(1,1)
[ydim,xdim,tdim]=size(fx);
sigma = zeros([ydim,xdim,tdim]);
for t = 1:tdim
    t
    if tt(t) == 0;
        
    else
        [dfxdx,dfxdy] = gradient(fx(:,:,t),dx,dy); 
        [dfydx,dfydy] = gradient(fy(:,:,t),dx,dy); 
        for i = 1:ydim
            for j = 1:xdim
                if ~isnan(dfxdx(i,j))&&~isnan(dfxdy(i,j))&&~isnan(dfydx(i,j))&&~isnan(dfydy(i,j))
                    gradF = [dfxdx(i,j),dfxdy(i,j); 
                             dfydx(i,j),dfydy(i,j)];
                    C = gradF'*gradF;
                    lambda = max(eig(C));
                    sigma(i,j,t) = 1/(2*tt(t))*log(lambda);
                else
                    sigma(i,j,t) = NaN;
                end
            end
        end
    end
end
time = tt;
save qge_sigma sigma time
