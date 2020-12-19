function MainPrgm
clc; clear all 
tic
clear t
clear y ys
clc;clf
global ysr y1r ysf y1f i j k  delta A1 B1 C1 D1 Ar1 Br1 Cr1 Dr1 Vd_on vin r
fs=5001; 
ts=1/fs; tmax=51*ts; delta=ts/503;
r=4;  % Desigred ordereduced
time=0:delta:tmax;
ysf=zeros(6,length(time));
ysr=zeros(6,length(time));
y1f=zeros(6,length(time));
y1r=zeros(r,length(time));
%% Full Model matrices
Vd_on = 0.7; vin=50;
RL=0.1; L=200e-6;   Ls=20e-9;  Rs=0.2; Cs=200e-12; Cd=100e-12; Rd=40e6; Rc=0.4; Lc=100e-12; C=42e-6;
Ro= 20;

A1 = [ (-RL - Ro)/L  0  Ro/L  0  Ro/L  -1/L ;
             0  0  0  0  1/C 0;
             Ro/Ls  0  -Ro/Ls  -1/Ls  -Ro/Ls  1/Ls;
             0  0  1/Cs  -1/(Cs*Rs)  0  0;
             Ro/Lc  -1/Lc  -Ro/Lc  0  (-Rc-Ro)/Lc  0;
             1/Cd  0  -1/Cd  0  0  -1/(Cd*Rd)];
B1 = [1/L  0 ;
    0  0;
    0  0 ;
    0  0;
    0  0;
    0 0/(Cd*Rd)];


C1 = eye(6);


D1 = zeros(6,2);

%% Reduced Model system matrices
sysFull = ss(A1,B1,C1,D1);

Wc = gram(sysFull,'c'); % Controllability Gramian
Wo = gram(sysFull,'o'); % Observability Gramian
%% Manually compute scaled balancing transformation
[Tu,D] = eig(Wc*Wo); % Tu are unscaled e-vecs

Atu = inv(Tu)*A1*Tu;
Btu = inv(Tu)*B1;
Ctu = C1*Tu;
Dtu = 0;
syst = ss(Atu,Btu,Ctu,Dtu);

Sigmac = gram(syst,'c'); % Diagonal Gramians
Sigmao = gram(syst,'o'); % (but not equal)
Sigmas = diag(Sigmac)./diag(Sigmao);

% Scaled balancing transformation
T = Tu*diag(Sigmas.^(1/4));

% Permute columns of T to order Sigma
Sigma = diag(Sigmac).^(1/2).*diag(Sigmao).^(1/2);
[sigsort,permind] = sort(Sigma,'descend');
T = T(:,permind); % Hierarchical

% Compute balanced system
At = inv(T)*A1*T;
Bt = inv(T)*B1;
Ct = C1*T;
Dt = 0;
sysBal = ss(At,Bt,Ct,Dt);

BWc = gram(sysBal,'c'); % Balanced Gramians
BWo = gram(sysBal,'o');

%%
sysBT = balred(sysFull,r);  % Balanced truncation
Ar1=sysBT.A;
Br1=sysBT.B;
Cr1=sysBT.C;
Dr1=sysBT.D;





%% Singular value Distribution of Full and Reduced systems
figure(1)
svdA1=svd(A1);
svdAr1=svd(Ar1);
subplot(211)
semilogy(svdA1,'b-*')
title('Singular values of distribution Full order sys.')
xlabel('i')
ylabel('\sigma_{i}')
grid on
subplot(212)
semilogy(svdAr1,'b-o')
grid on
title('Singular values of distribution Reduced order sys.')
xlabel('i')
ylabel('\sigma_{i}')

%%
i=1;
j=1;
k=1;
l=1;
tful=tic;
while k<=(length(time))
     modeful();
   k=k+1;
end
toc(tful)
%
tred=tic;
 while l<=(length(time))
     modereduced();
   l=l+1;
 end
 toc(tred)
%
 %-------------------FULL MODEL-------------------------
Yf=ysf';
aff=1:length(Yf(:,1));
%------------------------REDUCED----------------------------------
Yr=ysr';
arr=1:length(Yr(:,1));
 %% Comparision Plots of Reduced Model and Full Model
 %%%%%%%%%%%%%%%%%%%  Compare two Full and Reduced model dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 err=norm(Yf(:,1)-Yr(:,1));
 fprintf('Norm error of approximation is %f',err)
 figure(2)
  subplot(221)
    plot(aff, Yf(:,1),'b',arr,Yr(:,1),'r--') 
    hold on
    grid on
    title('Inductor Current')
   legend('Full','Reduced')
    subplot(222)
   plot(aff, Yf(:,2),'b',arr,Yr(:,2),'r--') 
    title('Diode Capacitor Voltage')
    legend('Full','Reduced')
    grid on
       subplot(223)
   plot(aff, Yf(:,3),'b',arr,Yr(:,3),'r--') 
    grid on
       title('Diode inductor current')
        legend('Full','Reduced')
       subplot(224)
   plot(aff, Yf(:,4),'b',arr,Yr(:,4),'r--') 
    grid on
       title('Capacitor Voltage')
       legend('Full','Reduced')
 
end
