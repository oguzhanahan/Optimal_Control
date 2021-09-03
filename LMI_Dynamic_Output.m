function [Ak,Bk,Ck,Dk] = LMI_Dynamic_Output(A,B1, B2, C1, C2, D11, D12, D21)


[n,n]=size(A)

[n,M]=size(B1)

[n,m]=size(B2)

[p2,M]=size(D21)

Xx=sdpvar(n)
Y=sdpvar(n)

Ac=sdpvar(n,n,'full')
Bc=sdpvar(n,p2,'full')
Cc=sdpvar(m,n,'full')
Dc=sdpvar(m,p2,'full')
gamma=sdpvar(1,1)


G11=Xx*A+A'*Xx+Bc*C2+C2'*Bc';
G12=Ac+A'+C2'*Dc'*B2'
G13=Xx*B1+Bc*D21
G14=C1'+C2'*Dc'*D12'

G21=G12';
G22=A*Y+Y*A'+B2*Cc+Cc'*B2';
G23=B1+B2*Dc*D21;
G24=Y*C1'+Cc'*D12';

G31=G13';
G32=G23';
G33=-gamma*eye(M);
G34=D11'+D21'*Dc'*D12';

G41=G14';
G42=G24';
G43=G34';
G44=-gamma*eye(6);


LMI1=[G11 G12 G13 G14
        G21 G22 G23 G24
        G31 G32 G33 G34
        G41 G42 G43 G44]<-.000001

LMI2=[Xx eye(n)
       eye(n) Y]>0.0000001

LMI3=trace(Xx)<1000
   Fset=[LMI1 LMI2 LMI3]
solution=optimize(Fset,gamma)   


Xx=value(Xx);
Y=value(Y);
Dc=value(Dc);
Cc=value(Cc);
Bc=value(Bc);
Ac=value(Ac);

Z=Y-Xx^-1;

Dk=Dc;
Ck=-Cc/Z+Dk*C2/Z;
Bk=Xx^-1*Bc-B2*Dk;
Ak=-Xx^-1*Ac*Z^-1+(A+Bk*C2)*Y*Z^-1-B2*Ck;



