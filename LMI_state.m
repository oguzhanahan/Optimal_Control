function K=LMI_state(A,B1,B2,C1,D11,D12)

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