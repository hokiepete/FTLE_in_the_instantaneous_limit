function dydt = dg_int(t,Y)
    %dydt = [2+tanh(Y(2)),0]';
    dydt = [2+tanh(Y(2)),-Y(2)-Y(2)^3]';
end

