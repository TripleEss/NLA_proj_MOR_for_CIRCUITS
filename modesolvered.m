function modesolvered(A,B,C,D,Voi)
%This is solver for reduced system ,
% Here the euler streatgry is used to solve the model in properformat
global y1r ysr i delta r
I=eye(r);
y1r(:,i+1)=(I-delta*A)\(y1r(1:end,i)+delta*B*Voi);
ysr(:,i+1)= C*y1r(:,i+1)+D*Voi;
end
