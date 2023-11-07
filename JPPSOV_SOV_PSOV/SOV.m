        clc;
        clear;
        close all;

        N=5000;
        sum_e = 0;
        time=100;            
        learning_rate_h=0.5;  %sov  h
for t=1:time
        t
        %x=rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        x=1*randn(1,N);x=(x-min(x))/(max(x)-min(x));  %w=0.5 h=0.3  x(kk).^2 (h过小收敛慢，h过大会有脉冲且稳态误差会增加)

        d=nonlinear1(x,N);
        d=Noise(20,d,N);
                    
        module_number=10;
        lrh=learning_rate_h;
        % 模块数
        M=module_number;
        % 每个模块的输出
        y=zeros(1,N);
        % 每个模块的输入个数
        N1=1+M+(M+1)*M/2;
        % 非线性部分的权重矩阵
        h=0.001*randn(1,N1);
        
        for i=M:1:N        %输入信号x的长度
      
            external_x=[x(i:-1:i-M+1)];

                % Volterra expansion  
                X_pre=external_x;
                L=M;
                X_post=X_pre;
                for ii=1:L
                   for jj=ii:L
                      X_post=[X_post X_pre(ii)*X_pre(jj)];
                   end
                end  
            X_post=[1 X_post].';  

            y(i)=h*X_post; 
            e(i)=d(i)-y(i);
            %非线性部分权重更新
            h=h+lrh*e(i)* X_post.'/((1+norm(X_post).^2));
                      
        end
         sum_e= sum_e + e.^2;
end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'c');
legend('SOV');


