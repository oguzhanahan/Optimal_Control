s=1
for I=1:6
    
subplot(3,2,s)
plot(t,x_unc.signals.values(:,I),'k--','LineWidth',2)
hold on
plot(t,x_cd.signals.values(:,I),'LineWidth',1)
hold on
plot(t,x_cc.signals.values(:,I),'r--','LineWidth',1)
met1=num2str(I);
met2='. Kat Yer Degisimi [m]';
legend('Kontrolsüz','Ayrik Kontrollü','Sürekli Kontrollü')
ylabel(strcat(met1,met2))
xlabel 'zaman [s]'
ylim([-0.5 1.5])
grid
s=s+1;
end



figure
plot(t,ud,'r--','LineWidth',4)
hold on
plot(t,uc,'k','LineWidth',4)
xlabel ('zaman [s]')
ylabel ('Kontrol Kuvveti [N]')
legend('Ayrik Kontrolcü','Sürekli Kontrolcü')
grid
