s=1
for I=1:6
 opts = bodeoptions('cstprefs');
 opts.PhaseVisible = 'off';
 opts.FreqUnits = 'Hz';
 opts.Xlim=[5 20];
 opts.Ylim=[-50 50];
    subplot(3,2,s)
sys1=ss(Ac,Bc1(:,1),Cc1(I,:),Dc11(I,1));
% sys2=ss(Ac+Bc2*Kc,Bc1(:,1),Cc1(I,:),Dc11(I,1));
sys3=ss(A+B1*Kd,E1(:,1),C(I,:),E2(I,1),T_s);
bodeplot(sys1,sys3,opts)
legend('Kontrolsüz','Ayrýk Kontrollü','Sürekli Kontrollü')
grid
met1=num2str(I);
met2='. Kat Yer Degisim Frekans Cevabi';
title(strcat(met1,met2));
xlabel 'Frekans';
ylabel 'Genlik';
s=s+1;
end




