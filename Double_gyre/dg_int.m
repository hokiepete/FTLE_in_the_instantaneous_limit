function dydt = dg_int(t,Y)
    A = 0.1;
    w = 0.2.*pi;
    e = 0.25;
    a = e.*sin(w.*t);
    b = 1-2.*e.*sin(w.*t);
    f = a.*Y(1).^2+b.*Y(1);
    dfdx = 2.*a.*Y(1)+b;    
    u =-pi.*A.*sin(pi.*f).*cos(Y(2).*pi);    
    v = pi.*A.*cos(pi.*f).*sin(Y(2).*pi).*dfdx;
    w = zeros(size(u));
    dydt = [u,v,w]';
end

