
%͹���ϵ���̶��벻�̶��ıȽ�
% load('PNFIR\data\const_PNFIR.mat','avg_e');
% const_e=avg_e;
% load('PNFIR\data\PNFIR0.01_0.6_0.6.mat','avg_e');
% e=avg_e;
% plot(1:5000, 10*log10(const_e), 'y',1:5000, 10*log10(e), 'r');
% legend('PNFIR(͹���ϵ��Ϊ1)','PNFIR')
%% 


%EFASTPNFIR��PNFIR�ıȽ�
% load('PNFIR\data\EFAST_PNFIR_L=1.mat','avg_e');
% EFAST_PNFIR=avg_e;
% load('PNFIR\data\PNFIR_L=1.mat','avg_e');
% PNFIR=avg_e;
% plot(1:5000, 10*log10(EFAST_PNFIR), 'y',1:5000, 10*log10(PNFIR), 'r');
% legend('EFAST_PNFIR','PNFIR')
%% 

%JPPSOBV��PNFIR������������бȽ�
% load ('JPPSOV_SOV_PSOV\data\JPPSOV0.2_0.5.mat','avg_ej');
% ej = avg_ej;
% load('PNFIR\data\rand\PNFIR0.01_0.6_0.6.mat','avg_e');
% e=avg_e;
% plot(1:5000, 10*log10(ej), 'b',1:5000, 10*log10(e), 'r');
% legend('JPPSOV','PNFIR')
%% 

%JPPSOBV��PNFIR����ɫ�������бȽ�
% load ('JPPSOV_SOV_PSOV\data\JPPSOV0.9_0.1.mat','avg_ej');
% ej = avg_ej;
% load('PNFIR\data\ARMA\PNFIR0.001_0.5_0.01.mat','avg_e');
% e=avg_e;
% plot(1:5000, 10*log10(ej), 'b',1:5000, 10*log10(e), 'r');
% legend('JPPSOV','PNFIR')


