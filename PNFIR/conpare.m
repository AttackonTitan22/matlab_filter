
%͹���ϵ���̶��벻�̶��ıȽ�
load('PNFIR\data\const_PNFIR.mat','avg_e');
const_e=avg_e;
load('PNFIR\data\PNFIR.mat','avg_e');
e=avg_e;
plot(1:5000, 10*log10(const_e), 'y',1:5000, 10*log10(e), 'r');
legend('PNFIR(͹���ϵ��Ϊ1)','PNFIR')