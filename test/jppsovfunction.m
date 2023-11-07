        clc;
        clear;
        close all;

        N=5000;
        sum_e = 0;
        time=100;     
        
        learning_rate_w=0.2;  %jppsov  w
        learning_rate_h=0.5;  %jppsov  h
for t=1:time
       t
       x=0.6*rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        %x=1*randn(1,N);x=(x-min(x))/(max(x)-min(x));  %w=0.5 h=0.3  x(kk).^2 (h过小收敛慢，h过大会有脉冲且稳态误差会增加)

        d=nonlinear1(x,N);   
        d=Noise(20,d,N);

        module_number=5;
        lrw=learning_rate_w;
        lrh=learning_rate_h;
        % Number of modules
        M=module_number;
        % Number of each module's input signed
        p=1;
        % Output of each module
        y=zeros(1,M);
        % Length of each module allover input
        N1=1+(p+1)+(p+1)*(p+2)/2;
        % Weight Vector matrix of nonlinear filter
        h=0.001*randn(1,N1);
        % Weight vector of linear section;       
        w=0.001*randn(1,M);
        % 每个模块的输入向量
        external_x=zeros(1,p+1);
        % 每个模块的输出延迟信号
        r=zeros(1,M);
        %wi(n)xi(n)
        
        y_final=zeros(N);
        
        e=JPPSOV_Model(N1,N,M,p,x,y,y_final,lrh,lrw,h,w,d);
        sum_e= sum_e + e.^2;
end



avg_ej = sum_e/time;
plot(1:5000, 10*log10(avg_ej), 'k');
legend('JPPSOV');

% save jppsov0.4.mat
