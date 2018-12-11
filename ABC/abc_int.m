function dydt = abc_int(t,y)
%ABC_FLOW Summary of this function goes here
  ABC_Amplitude=0.0;
  Ap = sqrt(3);
  Bp = sqrt(2);
  u = (Ap+ABC_Amplitude*sin(pi*t))*sin(y(3)) + cos(y(2));
  v = Bp*sin(y(1)) + (Ap+ABC_Amplitude*sin(pi*t))*cos(y(3));
  w = sin(y(2)) + Bp*cos(y(1));
  dydt = [u,v,w]';
end

