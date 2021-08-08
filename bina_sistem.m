clc
clear all
global A B1 B2


m1=450;
m2=345e3;
m3=345e3;
m4=345e3;


k1=18050;
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
     
     
 C=[c1+c2          -c2              0             0
       -c2                c2+c3      -c3            0 
         0                -c3              c3+c4    -c4 
         0                  0              -c4           c4]
     
     
     
  A=[zeros(4,4)         eye(4)
         -M^-1*K    -M^-1*C]
     
     
     
 B1=zeros(8,2)
 
 B1(5,1)=k1/m1;
 B1(5,2)=c1/m1;
 
 B2=[zeros(4,1);-1/m1;1/m2;0;0];
 
 C1=[eye(8),
          A(6,:)
          A(8,:)];

 
 D=[zeros(8,2);
        B1(6,:)
        B1(8,:)];
 
 [OV OD]=eig(A);
 
 OVV=eig(A)
%  plot(OVV,'o')
%  grid
%  xlabel 'Reel Eksen'
%  ylabel 'Ýmajiner Eksen'
% title('s düzlemi') 
%  
J(:,2)=imag(OVV);
 J(:,1)=real(OVV);
 
 