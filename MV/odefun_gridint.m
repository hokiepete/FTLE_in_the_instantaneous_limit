function dydt = odefun_gridint(t,Y,U,V)
    dydt = zeros(2,1);
    dydt(1) = U(Y(1),Y(2),t);%interp3(x,y,time,u,Y(1),Y(2),t,'spline');
    dydt(2) = V(Y(1),Y(2),t);%interp3(x,y,time,v,Y(1),Y(2),t,'spline');
end