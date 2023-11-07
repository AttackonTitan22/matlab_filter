        clc;
        clear;
        close all;

        N=5000;
        sum_e = 0;
        time=100;     
        
        learning_rate_w=0.03;  %  w越小收敛越慢
        learning_rate_h=0.001;  %  h>0.03 发散 


for t=1:time
       t
       %x=rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        x=1*randn(1,N);x=(x-min(x))/(max(x)-min(x));  %w=0.5 h=0.3  x(kk).^2 (h过小收敛慢，h过大会有脉冲且稳态误差会增加)
%          x=color_arma_signal_new(N);  %w=0.3 h=0.016  0.5*x(kk).^3
%         % x_a=0.36; 
%          x=(x-min(x))/(max(x)-min(x));
        d=zeros(1,N);
        
        % Nonlinear system identification
        for n=2:N
               d(n)=d(n-1)/(1+d(n-1).^2)+0.5*x(n).^3;
        end
   
      d=Noise(30,d,N);
                

        module_number=3;
        lrw=learning_rate_w;
        lrh=learning_rate_h;
        % 模块数
        M=module_number;
        % 外部输入个数
        L=2;
        % Output of each module
        z=zeros(N,M);
        y=zeros(N,M);
        % Weight Vector matrix of nonlinear filter
        h=0.001*randn(M,L+2);
        % Weight vector of linear section;       
        w=0.01*randn(1,M);
        % 每个模块的输入向量
        external_x=zeros(1,L+2);
        % 每个模块的输出延迟信号
        r=zeros(1,M);
        % 每个模块的反馈信号
        g=zeros(1,M);
        A=zeros(M,L+2);
        
        e=PNIIR1_Model(L,M,N,A,y,x,d,lrh,lrw,h,w);
        sum_e= sum_e + e.^2;
end


avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'k');
legend('PNIIR');

save ('PNIIR\data\PNIIR.mat','avg_e');



