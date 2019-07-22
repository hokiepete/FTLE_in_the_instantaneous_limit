function dydt = odefun_gridint(t,Y,U,V)
    dydt = zeros(2,1);
    dydt(1) = U(Y(1),Y(2),t);
    dydt(2) = V(Y(1),Y(2),t);
end
