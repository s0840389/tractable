close all
clear all

printirf=false;

pars

main_NK

main_YN

main_NKcap

main_YNcap

%% Figure IRF-RANK

figure(1)
clf

%output
subplot(3,2,1)

plot(IRF_NK(23,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YN(23,1:end-1)*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Output gap')

legend('NK','NK-YN','Location','southeast')

ylabel('PP')

%inflation
subplot(3,2,2)

plot(IRF_NK(11,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YN(11,1:end-1)*100,'r','LineWidth',1.8)



hline=refline(0,0);
hline.Color='black';
title('Inflation')

%consumption and investment
subplot(3,2,3)

plot(IRF_NK(19,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YN(19,1:end-1)*100,'r','LineWidth',1.8)

plot(IRF_NK(17,1:end-1)*100,'black--','LineWidth',1.8)
hold on

plot(IRF_YN(17,1:end-1)*100,'r--','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Consumption and Investment (dashed)')

ylabel('PP')

%Real wages
subplot(3,2,4)

plot(IRF_NK(6,2:end)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YN(6,2:end)*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Real wages')

%Labour share
subplot(3,2,5)

plot(IRF_NK(21,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YN(21,1:end-1)*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Labour share')
ylabel('PP')

%Productivity
subplot(3,2,6)

plot((IRF_NK(23,1:end-1)-IRF_NK(14,2:end))*100,'black','LineWidth',1.8)
hold on

plot((IRF_YN(23,1:end-1)-IRF_YN(14,2:end))*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Labour productivity')


%% Figure IRF-WCNK

figure(2)
clf

%output
subplot(3,2,1)

plot(IRF_NKcap(22,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YNcap(22,1:end-1)*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Output gap')

legend('NK','NK-YN','Location','southeast')

ylabel('PP')

%inflation
subplot(3,2,2)

plot(IRF_NKcap(10,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YNcap(10,1:end-1)*100,'r','LineWidth',1.8)



hline=refline(0,0);
hline.Color='black';
title('Inflation')

%Consumption: Worker & Entrepeur
subplot(3,2,3)

plot(IRF_NKcap(25,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YNcap(30,1:end-1)*100,'r','LineWidth',1.8)

plot(IRF_NKcap(24,1:end-1)*100,'black--','LineWidth',1.8)
hold on

plot(IRF_YNcap(29,1:end-1)*100,'r--','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Consumption: Worker and Entrepreneur (dashed)')

ylabel('PP')

%Investment
subplot(3,2,4)

plot(IRF_NKcap(16,2:end)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YNcap(16,2:end)*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Investment')

%Labour share
subplot(3,2,5)

plot(IRF_NKcap(20,1:end-1)*100,'black','LineWidth',1.8)
hold on

plot(IRF_YNcap(20,1:end-1)*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Labour share')
ylabel('PP')

%Productivity
subplot(3,2,6)

plot((IRF_NKcap(22,1:end-1)-IRF_NKcap(13,2:end))*100,'black','LineWidth',1.8)
hold on

plot((IRF_YNcap(22,1:end-1)-IRF_YNcap(13,2:end))*100,'r','LineWidth',1.8)

hline=refline(0,0);
hline.Color='black';
title('Labour productivity')



