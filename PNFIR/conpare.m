
%凸组合系数固定与不固定的比较
load('PNFIR\data\const_PNFIR.mat','avg_e');
const_e=avg_e;
load('PNFIR\data\PNFIR.mat','avg_e');
e=avg_e;
plot(1:5000, 10*log10(const_e), 'y',1:5000, 10*log10(e), 'r');
legend('PNFIR(凸组合系数为1)','PNFIR')