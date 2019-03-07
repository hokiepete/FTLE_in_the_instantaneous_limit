close all
clear all
clc
load flow_map_gridint_v2_long

dx = xx(1,2)- xx(1,1)
dy = yy(2,1)- yy(1,1)
[leny,lenx,tlen]=size(fx);
T = time - time(1);
sigma = zeros([leny,lenx,tlen]);
for t = 1:tlen
    t
    %need to skip t0 to avoid division by zero
    %set FTLE to zero instead
    if T(t) == 0
        continue
    else
        [dfxdx,dfxdy] = gradient(squeeze(fx(:,:,t)),dx,dy);
        [dfydx,dfydy] = gradient(squeeze(fy(:,:,t)),dx,dy);
        for i = 1:leny
            for j = 1:lenx
                if ~isnan(dfxdx(i,j))&&~isnan(dfxdy(i,j))&&~isnan(dfydx(i,j))&&~isnan(dfydy(i,j))
                    gradF = [dfxdx(i,j),dfxdy(i,j);
                             dfydx(i,j),dfydy(i,j)];
                    C = gradF'*gradF;
                    lambda = max(eig(C));
                    sigma(i,j,t) = 1/(2*abs(T(t)))*log(lambda);
                else
                    sigma(i,j,t) = nan;
                end
            end
        end
    end
end

save('FTLE_wrf_long.mat', 'sigma', 'T');
