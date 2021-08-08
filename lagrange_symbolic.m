clc
clear all

for I=1:11
    str.x{I}='x';
    str.d1x{I}='d1x';
    str.d2x{I}='d2x';
    str.m{I}='m';
    str.k{I}='k';
    str.c{I}='c';
    str.delta{I}='delta';
end
for I=1:8
    str.L{I}='L';
end
for I=1:5
    str.xy{I}='xy';
end
str.x=genvarname(str.x);
str.d1x=genvarname(str.d1x);
str.d2x=genvarname(str.d2x);
str.m=genvarname(str.m);
str.k=genvarname(str.k);
str.c=genvarname(str.c);
str.L=genvarname(str.L);
str.delta=genvarname(str.delta);
str.xy=genvarname(str.xy);
for I=1:11
    symbol.X(I)=sym(str.x{I});
    symbol.d1X(I)=sym(str.d1x{I});
    symbol.d2X(I)=sym(str.d2x{I});
    symbol.Mm(I)=sym(str.m{I});
    symbol.Kk(I)=sym(str.k{I});
    symbol.Cc(I)=sym(str.c{I});
    symbol.D(I)=sym(str.delta{I});
end
for I=1:8
    symbol.LL(I)=sym(str.L{I});
end
for I=1:5
    symbol.xy(I)=sym(str.xy{I});
end
%% Çekiciye ait gövdedeki süspansiyon baðlantýsý deplasmanlarý
symbol.xi=symbol.X(1)+symbol.LL(1)*symbol.X(8);
symbol.xj=symbol.X(1)-symbol.LL(2)*symbol.X(8);
symbol.xk=symbol.X(1)-symbol.LL(3)*symbol.X(8);
%% Çekici-Trayler arasý süspansiyon
symbol.xq1=symbol.X(1)-symbol.LL(4)*symbol.X(8);
symbol.xq2=symbol.X(2)+symbol.LL(5)*symbol.X(9);
%% Trayler arka teker takýmý
symbol.xp=symbol.X(2)-symbol.LL(6)*symbol.X(9);
%% trayler arka teker takýmý gövdedeki süspansiyon baðlantý deplasmanlarý
symbol.xm=symbol.xp+symbol.LL(7)*symbol.X(10);
symbol.xn=symbol.xp-symbol.LL(8)*symbol.X(10);


%% süspansiyon baðýl hareketi
%çekici-trayler arasý
symbol.D(1)=symbol.xq1-symbol.xq2;
% çekici ön teker 
symbol.D(2)=symbol.xi-symbol.X(3);
symbol.D(3)=symbol.X(3)-symbol.xy(1);
%çekici orta teker
symbol.D(4)=symbol.xj-symbol.X(4);
symbol.D(5)=symbol.X(4)-symbol.xy(2);
%çekici arka teker
symbol.D(6)=symbol.xk-symbol.X(5);
symbol.D(7)=symbol.X(5)-symbol.xy(3);
%trayler ön teker
symbol.D(8)=symbol.xm-symbol.X(6);
symbol.D(9)=symbol.X(6)-symbol.xy(4);
%trayler arka teker
symbol.D(10)=symbol.xn-symbol.X(7);
symbol.D(11)=symbol.X(7)-symbol.xy(5);

%% türevlenebilir çarpanlarýn gösterimi
for I=1:11
    symbol.dXd(I)=symbol.d2X(I)*symbol.d1X(I);
    symbol.Di(I)=symbol.D(I)^2;
    symbol.Ci(I)=symbol.D(I)^2;
end

symbol.Di=symbol.Di';

%% M, K ve C matrislerinin eldesi
for I=1:10
for J=1:10
    
symbol.M(I,J)=symbol.Mm(I)*diff(diff(symbol.dXd(J),symbol.d1X(J)),symbol.d2X(I));
symbol.K(I,J)=diff(0.5*symbol.Kk*diff(symbol.Di,symbol.X(I)),symbol.X(J));
symbol.C(I,J)=diff(0.5*symbol.Cc*diff(symbol.Di,symbol.X(I)),symbol.X(J));
end
symbol.M
symbol.K
symbol.C
end

for I=1:10
for J=1:5
    
symbol.B(I+10,J)=-diff(0.5*symbol.Kk*diff(symbol.Di,symbol.X(I)),symbol.xy(J))/symbol.Mm(I);
symbol.B(I+10,J+5)=-diff(0.5*symbol.Cc*diff(symbol.Di,symbol.X(I)),symbol.xy(J))/symbol.Mm(I);
end
symbol.B
end

%% Parametrelerin belirlenmesi
values.Mp=[9785/2 26000/2 270 520 520 270 270 18311/2 251900/2 3300 0];
values.Kp=[155800 300e3 847000 967430 2e6 967430 2e6 155800 2e6 155800 2e6];
values.Cp=[200e2 10e3 0 27627 0 27627 0 44506 0 44506 0];
values.Lp=[1.2 3.6 4.8 4.134 6.973 4 0.685 0.7];

symbol.Mv=symbol.M;
symbol.Kv=symbol.K;
symbol.Cv=symbol.C;
symbol.Bv=symbol.B;
for I=1:10
for J=1:10
for i=1:11
    
symbol.Mv(I,J)=subs(symbol.Mv(I,J),str.m(i),values.Mp(i));
symbol.Kv(I,J)=subs(symbol.Kv(I,J),str.k(i),values.Kp(i));
symbol.Cv(I,J)=subs(symbol.Cv(I,J),str.c(i),values.Cp(i));



end

for j=1:8
    
symbol.Kv(I,J)=subs(symbol.Kv(I,J),str.L(j),values.Lp(j));
symbol.Cv(I,J)=subs(symbol.Cv(I,J),str.L(j),values.Lp(j));

end

end

end

for I=1:20
    for J=1:10
        for i=1:11
symbol.Bv(I,J)=subs(symbol.Bv(I,J),str.k(i),values.Kp(i));
symbol.Bv(I,J)=subs(symbol.Bv(I,J),str.c(i),values.Cp(i));
symbol.Bv(I,J)=subs(symbol.Bv(I,J),str.m(i),values.Mp(i));
        end
end
end
values.Mv=double(symbol.Mv);
values.Kv=double(symbol.Kv);
values.Cv=double(symbol.Cv);
values.Bv=double(symbol.Bv);

values.Av=[zeros(size(values.Mv)) eye(size(values.Mv))
    -values.Kv*values.Mv^-1      -values.Cv*values.Mv^-1];

% [values.EV,values.ED]=eig(values.Av);
% 
% plot(diag(values.ED),'o');

values.Cv=[1 zeros(1,9)];
values.Dv=zeros(1,10);

clear i, clear j, clear I, clear J

%% insan modeli
  ms1=15;  ms2=1+7.8;  ms3=43.4;    
    ks1=31000;  ks2=18000;   ks3=44130;
    cs1=830; cs2=200;   cs3=1485;
    
    
    As=[    0                        1                   0               0                              0                          0
            -ks1/ms1      -(cs1+cs2)/ms1   ks2/ms1       cs2/ms1                    0                          0
               0                       -1                  0               1                               0                          0
               0                     cs2/ms2          -ks2/ms2     -(cs2+cs3)/ms2       ks3/ms2               cs3/ms2
               0                        0                  0              -1                               0                          1
               0                        0                  0               cs3/ms3                   -ks3/ms3             -cs3/ms3];
           
  %% x_dot=Ax+B1w+B2u         
           
Ac=[values.Av                                           zeros(length(values.Av),length(As))
       zeros(length(As),length(values.Av))                         As                   ];               
r=length(values.Av );
Ac(r+1,r/2+1)=-1;

Bs=zeros(6,1);
Bs(2,1)=-1/ms1;

Bc1=[ values.Bv
          zeros(6,10)]
      
 Bc2=[zeros(length(values.Av),1)
           Bs]

 %% z=C1x+D11w+D12u
 Cc1=[zeros(6,length(Ac)-6) eye(6)]
 
 Dc11=zeros(6,10)
 
 Dc12=zeros(6,1)
 
 %% y=C2x+D21w
 Cc2=eye(length(Ac))
 
 Dc21=zeros(size(Cc2,1),10)