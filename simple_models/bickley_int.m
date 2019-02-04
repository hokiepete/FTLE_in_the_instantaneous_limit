function dydt = dg_int(t,Y)
    U0 = 62.66;%m s^-1
    U0 = 24*3.6*U0; %km daY(2)^-1
    L = 1770;%km
    A2 = 0.1;
    A3 = 0.3;
    c2 = 0.205*U0;
    c3 = 0.461*U0;
    re = 6371;%km
    k = @(n) 2*n/re;
    dpdx = -k(3).*A3.*L.*U0.*((sech(Y(2)./L)).^2).*sin(k(3).*(Y(1)-c3.*t)) - k(2).*A2.*L.*U0.*((sech(Y(2)./L)).^2).*sin(k(2).*(Y(1)-c2.*t));
    dpdy = - U0.*((sech(Y(2)./L)).^2) - 2.*A3.*U0.*tanh(Y(2)./L).*((sech(Y(2)./L)).^2).*cos(k(3).*(Y(1)-c3.*t)) - 2.*A2.*U0.*tanh(Y(2)./L).*((sech(Y(2)./L)).^2).*cos(k(2).*(Y(1)-c2.*t));
    dydt = [-dpdy,dpdx]';
end

