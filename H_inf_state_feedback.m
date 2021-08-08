clc
clear all

Pen=5e-6;

m1=350;
m2=10;
k2=200000;
k1=15000;
c1=335;


global A B1 B2

A=[0 0                                                1 0;
        0 0                                                0 1;
        -k1/m1 k1/m1                       -c1/m1 c1/m1;
        k1/m2 -k1/m2-k2/m2             c1/m2 -c1/m2];

B1=[0;0;0;k2/m2];

B2=[0;0;1/m1;-1/m2];

C1=[1 0 0 0;
         0 1 0 0
         0 0 0 0
         0 0 0 0];



D11=zeros(size(C1,1),size(B1,2));
 
D11(1,:)=Pen;

D12=zeros(size(C1,1),size(B2,2)) ;

D12(2,:)=Pen;

[n m]=size(B2);
p=size(B1,2);
c=size(C1,1);


X=sdpvar(n,n);
W=sdpvar(m,n,'full');
gamma=sdpvar(1,1);

LMI1=[A*X+X*A'+B2*W+W'*B2'           B1            X*C1'+W'*D12';
             B1'               -gamma*eye(p)          D11';
             C1*X+D12*W              D11           -gamma*eye(c)]<0;
         
LMI2=X>0;

Fset=[LMI1,LMI2];

solution=optimize(Fset,gamma);

W=value(W);
X=value(X);
gammac=value(gamma)
K=W*inv(X);

olp=eig(A)
clp=eig(A+B2*K)

figure;

plot(olp,'o')
hold on
grid on
plot(clp,'or')

sim('simulate_H_inf_state')

figure
Time_Domain_Plots
figure
Frequency_Domain_Plots