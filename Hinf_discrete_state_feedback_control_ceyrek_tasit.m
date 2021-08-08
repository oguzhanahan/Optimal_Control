clc
clear all

T_s=1e-2;

global A B1 E1 

m1=350;
m2=10;
k2=500000;
k1=10000;
c1=9000;

Ac=[0 0 1 0;
        0 0 0 1;
        -k1/m1 k1/m1 -c1/m1 c1/m1;
        k1/m2 -k1/m2-k2/m2 c1/m2 -c1/m2];
Bc1=[0;0;0;k2/m2];
Bc2=[0;0;1/m1;-1/m2];
Cc1=[1 0 0 0;0 1 0 0];

  
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
% 
% 
% C=[1  0 0 
%         0   1 0];
% 
% B2=[0 
%         0];
% 
% E2=[0 
%         0];



% Ac=[-1.01887           0.90506           -0.00215
%               0.82225      -1.07741          -0.17555
%                 0                   0                   -20.2];
            
[ar ac]=size(Ac)

% Bc1=[0
%          0
%          20.2];
     
[b1r b1c]=size(Bc1)

%      Bc2=[1
%                0
%                0];



% Cc1=[1  0 0 
%         0   1 0];

[c1r c1c]=size(Cc1)

Dc11=[0 
        0];

Dc12=[0 
            1];


Cc2=eye(c1c)

Dc21=zeros(c1c,b1c)
    
Dc22=zeros(c1c,b1c)


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

B2=d(1:size(Cc1,1),1:size(Bc1,2))

E2=d(1:size(Cc1,1),size(Bc1,2)+1:size(b,2))


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
           G41   G42      G43         G44]<-0.00001
  
  
  LMI2=Q>0;
  
  
  
  Fset=[LMI1   LMI2]
  
  solution=optimize(Fset,gamma)
  
  

Q=value(Q);
Y=value(Y);
  
K=Y*Q^-1

eig(Q)

sys_d_c=ss(A+B1*K,E1,[C+B2*K;eye(c1c)+zeros(c1c,b1c)*K],[E2;zeros(c1c,b1c)],T_s);

eig(A+B1*K)

ltiview(sys_d_c)
