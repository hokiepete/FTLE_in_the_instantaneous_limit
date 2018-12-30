function dydt = abc_flow(x,y,z,t)
%ABC_FLOW Summary of this function goes here
  ABC_Amplitude=0.0;
  Ap = sqrt(3);
  Bp = sqrt(2);
  u = (Ap+ABC_Amplitude*sin(pi*t))*sin(z) + cos(y);
  v = Bp*sin(x) + (Ap+ABC_Amplitude*sin(pi*t))*cos(z);
  w = sin(y) + Bp*cos(x);
  dydt = [u,v,w]';
end

