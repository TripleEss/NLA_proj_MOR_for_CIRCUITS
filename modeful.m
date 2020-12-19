function ggf =modeful()
% This is the sub program which is called from main program and itself it
% calls the modesolvelful, calls solver for Full order  system
global j  A1 B1 C1 D1 Vd_on vin
u1=[vin;Vd_on];
modesolvefull(A1,B1,C1, D1,u1)
j=j+1;
end