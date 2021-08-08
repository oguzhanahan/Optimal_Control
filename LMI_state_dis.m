function K=LMI_state_dis(A,B1,E1,C,B2,E2)

[n,n]=size(A);

[n,m]=size(B1);

[n,r]=size(E1);

[q,n]=size(C);
%% LMI 


Q=sdpvar(n,n);
Y=sdpvar(m,n,'full')
gamma=sdpvar(1,1)


G11=-Q;                 G12=zeros(n,r);                           G13=Q*A'+Y'*B1';                 G14=Q*C'+Y'*B2';

G21=G12';              G22=-gamma*eye(r);                    G23=E1';                                  G24=E2';


G31=G13';               G32=G23';                                  G33=-Q;                                   G34=zeros(n,q);

G41=G14';               G42=G24';                                  G43=G34';                                G44=-gamma*eye(q);


 
 
 
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
K=Y*Q^-1