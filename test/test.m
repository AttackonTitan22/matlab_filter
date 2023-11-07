        clc;
        clear;
        close all;

        N=5000;
        sum_e = 0;
        time=10;     
        
        learning_rate_w=0.1;  %  w
        learning_rate_h=0.1;  %  h

        
        TT=1;
        Ta=0;
for t=1:time
       t
       x=rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        %x=0.6*randn(1,N);x=(x-min(x))/(max(x)-min(x));  %w=0.5 h=0.3  x(kk).^2 (h过小收敛慢，h过大会有脉冲且稳态误差会增加)
%          x=color_arma_signal_new(N);  %w=0.3 h=0.016  0.5*x(kk).^3
%         % x_a=0.36; 
%          x=(x-min(x))/(max(x)-min(x));
        d=zeros(1,N);
        % Nonlinear system identification
        for n=2:N
               d(n)=d(n-1)/(1+d(n-1).^2)+x(n).^3;
        end
   
        snr=30;
        info_power=(norm(d)^2)/N;  %信号的能量                  
        snr_lin=10^(snr/10);                         
        noise_power=info_power/snr_lin;   %噪声的方差
        d=d+sqrt(noise_power)*randn(1,N);
        %    Parameters Lists:
        %                       train_x:  trainning sample;
        %                  train_number:  trainning number;
        %                 module_number:  module's number;
        %                  input_number:  input signal length 's number;
        %                      desire_d:  desire of output signal through the unknown system; 
        %                      lt:    parameter list;
        %   Output results variables:
        %                    e_d:  Error of output; 
        %                y_final:  final output  through the processing to the data;
        %                       
        train_number=N;
        module_number=5;
        lrw=learning_rate_w;
        lrh=learning_rate_h;
        % Number of modules
        M=module_number;
        % Length of external input
        L=2;
        % Length of trainning sampling.
        n=train_number;
        % Output of each module
        z=zeros(N,M);
        y=zeros(N,M);
        % Weight Vector matrix of nonlinear filter
        h=0.001*randn(M,L+1);
        % Weight vector of linear section;       
        w=0.001*randn(1,M);
        % 每个模块的输入向量
        external_x=zeros(1,L+1);
        % 每个模块的输出延迟信号
        r=zeros(1,M);
        
        for i=L+M-1:1:N        %输入信号x的长度
            A=zeros(M,L+1);
            for j=M:-1:1    %模块的长度
                if j==M
                    r(j)=y(i-1,1);
                else
                    r(j)=y(i-1,j+1);
                end
                external_x=[x(i-j+1:-1:i-(j+L-2)) r(j)]; 
                z(i,j)=h(j,:)*external_x.';
                y(i,j)=1/(1+exp(-z(i,j)));
               % A(j,:)=(TT*w(j)*(exp(-h(j,:)*external_x.')/(exp(-h(j,:)*external_x.') + 1)^2)+(1-TT))*external_x;
                A(j,:)=w(j)*external_x*(exp(-h(j,:)*external_x.')/(exp(-h(j,:)*external_x.') + 1)^2);
            end
            y_final(i)=w*y(i,:).';
            e(i)=d(i-1)-y_final(i); 
                %非线性部分权重更新
            h(j,:)=h(j,:)+lrh*e(i)* A(j,:)/((norm(A(j,:)).^2));
                 %线性部分权重更新
            w=w+lrw*e(i)*y(i,:)/(1+norm(y(i,:)).^2); 
               
        end
        sum_e= sum_e + e.^2;
end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'k');
legend('PNFIR(凸组合系数为1)');
% save const_PNFIR.mat

