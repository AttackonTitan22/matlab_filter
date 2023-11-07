%JPPSOV与PNFIR对比图
load PNFIR\data\jppsov0.4.mat
je=avg_ej;
load PNFIR\data\PNFIR.mat
e=avg_e;
plot(1:5000, 10*log10(je), 'b',1:5000, 10*log10(e), 'r');
legend('JPPSOV)','PNFIR')