function modesolvefull(A,B,C,D,Voi)
%This is solver for Full order system
% Here the euler streatgry is used to solve the model in properformat
global y1f ysf j delta
I=eye(4);
y1f(:,j+1)=(I-delta*A)\(y1f(1:end,j)+delta*B*Voi);
ysf(:,j+1)= C*y1f(:,j+1)+D*Voi;
end