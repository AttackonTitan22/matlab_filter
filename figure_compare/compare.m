%jppsov与sov与psov对比图
load ('JPPSOV_SOV_PSOV\data\sov.mat','avg_e');
load ('JPPSOV_SOV_PSOV\data\jppsov.mat','avg_ej');
load ('JPPSOV_SOV_PSOV\data\psov.mat','avg_ep');
plot(1:5000, 10*log10(avg_e),'c');
hold on
plot( 1:5000, 10*log10(avg_ej), 'k');
hold on
plot(1:5000, 10*log10(avg_ep), 'm');
hold off
legend('SOV','JPPSOV','PSOV');