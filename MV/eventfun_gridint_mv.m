function [value,isterminal,direction]=eventfun_gridint_mv(t,Y,U,V)
    value=[isnan(Y(1))-1,isnan(Y(2))-1];
    isterminal=[1,1];
    direction=[0,0];
end
%{

function [value,isterminal,direction]=eventfun_gridint_mv(t,Y,U,V)
    if isnan(Y(1))
        val_1 = 1e-14;
    else
        val_1 = -1e14;
    end
    if isnan(Y(2))
        val_2 = 1e-14;
    else
        val_2 = -1e-14;
    end
    value=[val_1,val_2]
    isterminal=[1,1];
    direction=[0,0];
end
%}