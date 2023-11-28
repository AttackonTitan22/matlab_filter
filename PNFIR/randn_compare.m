
%测试版用于寻找最佳对比图
%EFAST_PNFIR和PNFIR对比图  
load ('PNFIR\data\randn\PNFIR_test.mat','avg_e');
PNFIR=avg_e;
load ('PNFIR\data\randn\EFASTPNFIR_test.mat','avg_e');
EFAST_PNFIR=avg_e;
plot(1:5000, 10*log10(PNFIR), 'g',1:5000, 10*log10(EFAST_PNFIR), 'r');
legend('PNFIR','EFAST  PNFIR')



%EFAST_PNFIR和PNFIR对比图  w=0.1 h=0.8 a=0.6  λ=1  Ta=0
% load ('PNFIR\data\randn\PNFIR_0.1_0.8_0.6.mat','avg_e');
% PNFIR=avg_e;
% load ('PNFIR\data\randn\EFASTPNFIR_0.1_0.8_0.6.mat','avg_e');
% EFAST_PNFIR=avg_e;
% plot(1:5000, 10*log10(PNFIR), 'g',1:5000, 10*log10(EFAST_PNFIR), 'r');
% legend('PNFIR','EFAST  PNFIR')

%EFAST_PNFIR和PNFIR对比图  w=1 h=1 a=0.1  λ=1  Ta=0
% load ('PNFIR\data\randn\PNFIR_1_1_0.1.mat','avg_e');
% PNFIR=avg_e;
% load ('PNFIR\data\randn\EFASTPNFIR_1_1_0.1.mat','avg_e');
% EFAST_PNFIR=avg_e;
% plot(1:5000, 10*log10(PNFIR), 'g',1:5000, 10*log10(EFAST_PNFIR), 'r');
% legend('PNFIR','EFAST  PNFIR')