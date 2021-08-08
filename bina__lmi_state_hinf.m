clc
clear all
% lagrange_symbolic   
elcent_v
elcent_x
s=0;
T_s=0.01
global A B1 E1 
 

m1=862.85;
m2=862.85;
m3=862.85;
m4=862.85;
m5=862.85;
m6=803.98;


k1=1.26e6;
k2=1.23e6;
k3=1.23e6;
k4=1.23e6;
k5=1.23e6;
k6=1.23e6;



M=[m1 0     0     0     0    0
       0   m2   0     0     0    0
       0    0    m3   0     0     0
       0    0    0    m4    0    0
       0    0    0     0    m5  0
       0    0    0    0    0      m6]
   
   
K=[k1+k2          -k2              0              0            0         0
       -k2                k2+k3      -k3            0            0         0
         0                -k3              k3+k4    -k4          0         0
         0                  0              -k4            k4+k5  -k5       0
         0                  0                0            -k5       k5+k6   -k6
         0                  0                0              0        -k6          k6      ]
     
     ksi1=0.2/100;
     ksi2=0.09/100;
     
     w1=1.05*2*pi;
     w3=6.12*2*pi;
     
     
     RP=(2*w1*w3/(w3^2-w1^2))*[w3    -w1; -w3^-1  w1^-1]*[ksi1;ksi2]

% RP(1,1)=2*0.09*w1*w3/(w1+w3)
% RP(2,1)=2*0.09/(w1+w3)
     
Cd=RP(1,1)*M+RP(2,1)*K


Ac=[zeros(length(K)) eye(length(K));
        -M^-1*K             -M^-1*Cd]

 Bc1=zeros(size(Ac,1),2)
 
 Bc1(size(M,1)+1,1)=k1/m1;
 Bc1(size(M,1)+1,2)=(Cd(1,1)-Cd(1,2))/m1;
 
 Bc2=[zeros(size(M,1),1);-1/m1;1/m2;zeros(size(M,1)-2,1)];
 
 Cc1=[eye(size(Ac,1))];


[b1r b1c]=size(Bc1)
[c1r c1c]=size(Cc1)

Pen=1e-4
Dc11=zeros(size(Cc1,1),size(Bc1,2))
Dc11(1,1)=Pen;
Dc11(2,1)=Pen;
Dc11(3,1)=Pen;
Dc11(4,1)=Pen;
Dc11(5,1)=Pen;
Dc11(6,1)=Pen;
pen=1e-4
Dc12=zeros(size(Cc1,1),size(Bc2,2))
Dc12(1,1)=pen;
Dc12(2,1)=pen;
Dc12(3,1)=pen;
Dc12(4,1)=pen;
Dc12(5,1)=pen;
Dc12(6,1)=pen;
% Dc12(7,1)=pen;
% Dc12(8,1)=pen;
% Dc12(9,1)=pen;
% Dc12(10,1)=pen;
% Dc12(11,1)=pen;
% Dc12(12,1)=pen;

Cc2=eye(size(Ac,1))

Dc21=zeros(size(Cc2,1),size(Bc1,2))
Dc22=zeros(size(Cc2,1),size(Bc2,2))

Kc=LMI_state(Ac,Bc1,Bc2,Cc1,Dc11,Dc12)

[A,B1,E1,C,B2,E2]=dis_state_space(Ac,Bc1,Bc2,Cc1,Dc11,Dc12,Cc2,Dc21,Dc22,T_s)

Kd=LMI_state_dis(A,B1,E1,C,B2,E2)


sim('simulate_H_inf_state.slx')


% figure
% Time_Domain_Plots
