
%测试版用于寻找最佳对比图
%EFAST_PNFIR和PNFIR对比图  
load ('PNFIR\data\ARMA\PNFIR_test.mat','avg_e');
PNFIR=avg_e;
load ('PNFIR\data\ARMA\EFASTPNFIR_test.mat','avg_e');
EFAST_PNFIR=avg_e;
plot(1:5000, 10*log10(PNFIR), 'g',1:5000, 10*log10(EFAST_PNFIR), 'r');
legend('PNFIR','EFAST  PNFIR')

%EFAST_PNFIR和PNFIR对比图  w=0.5 h=0.9 a=0.5  λ=1   Ta=0
% load('PNFIR\data\ARMA\PNFIR_0.5_0.9_0.5.mat','avg_e');
% PNFIR=avg_e;
% load('PNFIR\data\ARMA\EFASTPNFIR_0.5_0.9_0.5.mat','avg_e');
% EFAST_PNFIR=avg_e;
% plot(1:5000, 10*log10(PNFIR), 'g',1:5000, 10*log10(EFAST_PNFIR), 'r');
% legend('PNFIR','EFAST  PNFIR')