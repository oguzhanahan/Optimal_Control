clc
clear all
% lagrange_symbolic   
T_s=0.01;
global A B1 E1 
 

m1=450;
m2=345e3;
m3=345e3;
m4=345e3;


k1=1805;
k2=340000;
k3=326000;
k4=280000;


c1=26170;
c2=490;
c3=467;
c4=410;


M=[m1 0     0     0
       0   m2   0     0
       0    0    m3   0
       0    0    0    m4]
   
   
K=[k1+k2          -k2              0              0
       -k2                k2+k3      -k3            0 
         0                -k3              k3+k4    -k4 
         0                  0              -k4            k4]
     
     
 Cd=[c1+c2          -c2              0             0
       -c2                c2+c3      -c3            0 
         0                -c3              c3+c4    -c4 
         0                  0              -c4           c4]

Ac=[zeros(4) eye(4);
        -M^-1*K -M^-1*Cd]

 Bc1=zeros(8,2)
 
 Bc1(5,1)=k1/m1;
 Bc1(5,2)=c1/m1;
 
 Bc2=[zeros(4,1);-1/m1;1/m2;0;0];
 
 Cc1=[eye(8)];


[b1r b1c]=size(Bc1)
[c1r c1c]=size(Cc1)

Pen=1e-3
Dc11=zeros(size(Cc1,1),size(Bc1,2))
Dc11(1,1)=Pen;
Dc12=zeros(size(Cc1,1),size(Bc2,2))

Cc2=eye(size(Ac,1))

Dc21=zeros(size(Cc2,1),size(Bc1,2))
Dc22=zeros(size(Cc2,1),size(Bc2,2))

Kc=LMI_state(Ac,Bc1,Bc2,Cc1,Dc11,Dc12)

[A,B1,E1,C,B2,E2]=dis_state_space(Ac,Bc1,Bc2,Cc1,Dc11,Dc12,Cc2,Dc21,Dc22,T_s)

Kd=LMI_state_dis(A,B1,E1,C,B2,E2)


sim('simulate_H_inf_state.slx')
figure
Time_Domain_Plots
