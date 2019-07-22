function [value,isterminal,direction]=eventfun_gridint_mv(t,Y,U,V)
    value=[isnan(Y(1))-1,isnan(Y(2))-1];
    isterminal=[1,1];
    direction=[0,0];
end