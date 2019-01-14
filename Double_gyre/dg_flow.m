function [u,v,w] = abc_flow(x,y,z,t)
    A = 0.1;
    w = 0.2.*pi;
    e = 0.25;
    a = e.*sin(w.*t);
    b = 1-2.*e.*sin(w.*t);
    f = a.*x.^2+b.*x;
    dfdx = 2.*a.*x+b;    
    u =-pi.*A.*sin(pi.*f).*cos(y.*pi);    
    v = pi.*A.*cos(pi.*f).*sin(y.*pi).*dfdx;
    w = zeros(size(u));
end

