function ggr =modereduced()
% This is the sub program which is called from main program and itself it
% calls the modesolvered, calls solver for reduced system
global i  Ar1 Br1 Cr1 Dr1 
Vd_on = 0.7;
vin=50; 
u1=[vin;Vd_on];
modesolvered(Ar1,Br1,Cr1, Dr1,u1)
i=i+1;
end