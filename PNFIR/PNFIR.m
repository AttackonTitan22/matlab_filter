        clc;
        clear;
        close all;

        N=5000;
        sum_e = 0;
        time=100;     
        
        learning_rate_w=0.01;  %  w
        learning_rate_h=0.6;  %  h
        learning_rate_a=0.6;  %  a
        
        TT=1;
        Ta=0;
for t=1:time
        t
        x=0.6*rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        %x=randn(1,N);x=0.4*(x-min(x))/(max(x)-min(x)); %w=0.5 h=0.3  x(kk).^2 (h��С��������h���������������̬��������)
        
        d=nonlinear1(x,N);
        d=Noise(30,d,N);
                  

        module_number=5;
        lrw=learning_rate_w;
        lrh=learning_rate_h;
        lra=learning_rate_a;

        % Number of modules
        M=module_number;
        % Length of external input
        L=2;
        % Output of each module
        z=zeros(N,M);
        y=zeros(N,M);
        % Weight Vector matrix of nonlinear filter
        % h=0.001*randn(M,L+1);
        h=0.001*randn(M,L+1);
        % Weight vector of linear section;       
        w=0.001*randn(1,M);
        % ÿ��ģ�����������
        X=zeros(1,L+1);
        % ÿ��ģ�������ӳ��ź�
        r=zeros(1,M);
        A=zeros(M,L+1);

        %����PNFIR����������
        e=PNFIR_Model(L,M,N,A,y,x,d,lrh,lrw,lra,h,w,TT,Ta);      
        sum_e= sum_e + e.^2;
end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'k');
legend('PNFIR');

save ('PNFIR\data\PNFIR.mat','avg_e');