clc
clear all
% lagrange_symbolic   
T_s=0.0001;
global A B1 E1 
 

m1=350;
m2=10;
k2=200000;
k1=15000;
c1=335;

K=[k1 -k1;-k1 k1+k2];

M=[m1 0; 0 m2];

Cd=[c1 -c1;-c1 c1]

Ac=[zeros(2) eye(2);
        -M^-1*K -M^-1*Cd]

Bc1=[0;0;0;k2/m2];

Bc2=[0;0;1/m1;-1/m2];

Cc1=[1 0 0 0;
         0 1 0 0
         0 0 0 0
         0 0 0 0];

%    x_k+1=A*x_k   +      B1*u_k    +    E1*w_k
% A=[       9.0649*10^-1      8.1601 * 10^-2      -5.0128 * 10^-4
%              7.4135 * 10^-2     9.0121 * 10^-1     -7.0423 * 10^-3
%               0                          0                              1.3266 * 10^-1];
%       
% 
%      B1=[  -4.4486*10^-4
%               -9.6285*10^-3
%                8.6734*10^-1];
%        
% E1=[9.5189 *10^-2
%         3.8373 *10^-3
%          0                    ]  ;
% z_k=C*x_k    +         B2*u_k       +      E2*w_k
% Ac=[-1.01887           0.90506           -0.00215
%               0.82225      -1.07741          -0.17555
%                 0                   0                   -20.2];
%             
% [ar ac]=size(Ac)
% 
% Bc1=[0
%          0
%          20.2];
%      
% 
%      Bc2=[1
%                0
%                0];
% 
% 
% 
% Cc1=[1  0 0 
%         0   1 0];

[b1r b1c]=size(Bc1)

[c1r c1c]=size(Cc1)

Pen=0
Dc11=zeros(size(Cc1,1),size(Bc1,2))
Dc11(1,1)=Pen;
Dc12=zeros(size(Cc1,1),size(Bc2,2))
Dc12(2,1)=Pen;

Cc2=eye(size(Ac,1))

Dc21=zeros(size(Cc2,1),size(Bc1,2))


Dc22=zeros(size(Cc2,1),size(Bc2,2))


AA=Ac;

BB=[Bc1     Bc2]

CC=[Cc1;Cc2]

DD=[Dc11    Dc12
         Dc21    Dc22]
     
     
sys=ss(AA,BB,CC,DD)

sys_d=c2d(sys,T_s,'zoh')

[a,b,c,d] = ssdata(sys_d)

A=a

B1=b(:,1:size(Bc1,2))
     
E1=b(:,size(Bc1,2)+1:size(b,2))

C=c(1:size(Cc1,1),:)

B2=d(1:size(C,1),1:size(Bc1,2))

E2=d(1:size(C,1),size(Bc1,2)+1:size(b,2))


[n,n]=size(A);

[n,m]=size(B1);

[n,r]=size(E1);

[q,n]=size(C);
%% LMI 


Q=sdpvar(n,n);
Y=sdpvar(m,n,'full')
gamma=sdpvar(1,1)


G11=-Q;                 G12=zeros(n,r);                           G13=Q*A'+Y'*B1';                 G14=Q*C'+Y'*B2';

G21=G12';              G22=-gamma*eye(r);            G23=E1';                                  G24=E2';


G31=G13';               G32=G23';                               G33=-Q;                                   G34=zeros(n,q);

G41=G14';               G42=G24';                               G43=G34';                                G44=-gamma*eye(q);


 
 
 
LMI1=[G11   G12       G13        G14
           G21   G22       G23        G24
           G31   G32       G33        G34
           G41   G42      G43         G44]<0
  
  
  LMI2=Q>0;
  
  
  
  Fset=[LMI1   LMI2]
  
  solution=optimize(Fset,gamma)
  
  

Q=value(Q);
Y=value(Y);
gamma=value(gamma)  
Kd=Y*Q^-1

sim('simulate_H_inf_state.slx')
figure
Time_Domain_Plots

