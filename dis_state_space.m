function [A B1 E1 C B2 E2]=dis_state_space(Ac,Bc1,Bc2,Cc1,Dc11,Dc12,Cc2,Dc21,Dc22,T_s)

AA=Ac;

BB=[Bc1     Bc2]

CC=[Cc1;Cc2]

DD=[Dc11    Dc12
         Dc21    Dc22]
     
     
sys=ss(AA,BB,CC,DD)

sys_d=c2d(sys,T_s,'zoh')

[a,b,c,d] = ssdata(sys_d)

A=a %x

E1=b(:,1:size(Bc1,2))  %w
     
B1=b(:,size(Bc1,2)+1:size(b,2)) %u

C=c(1:size(Cc1,1),:) %x

E2=d(1:size(C,1),1:size(Bc1,2)) %w

B2=d(1:size(C,1),size(Bc1,2)+1:size(b,2)) %u