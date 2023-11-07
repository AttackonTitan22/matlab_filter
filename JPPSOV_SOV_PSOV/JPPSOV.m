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
%          x=color_arma_signal_new(N);  %w=0.3 h=0.016  0.5*x(kk).^3
%         % x_a=0.36; 
%          x=(x-min(x))/(max(x)-min(x));
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
        external_s=zeros(1,p+1);
        % 每个模块的输出延迟信号
        r=zeros(1,M);
        %wi(n)xi(n)
        
        y_final=zeros(N);
        e=zeros(N);
        
        for i=N1:1:N        %输入信号x的长度
            temp_w_sum=zeros(1,N1);
            for j=M:-1:1    %模块的长度
                if j==M
                    r(j)=y(j);
                else
                    r(j)=y(j+1);
                end
                external_x=[x(i-j:-1:i-(j+p-1)) r(j)];
                    % Volterra expansion  
                    X_pre=external_x;
                    L=p+1;
                    X_post=X_pre;
                    for ii=1:L
                       for jj=ii:L
                          X_post=[X_post X_pre(ii)*X_pre(jj)];
                       end
                    end  
                X_post=[1 X_post].';  
                y(j)=h*X_post;
                temp_w_sum=temp_w_sum+w(j)*X_post.';    

            end
            y_final(i)=w*y.';
            e(i)=d(i-1)-y_final(i);
            %线性部分权重更新
            w=w+lrw*e(i)*y/(1+norm(y).^2);  
            %非线性部分权重更新
            h=h+lrh*e(i)* temp_w_sum/((1+norm(temp_w_sum).^2));  
        end
        sum_e= sum_e + e.^2;
end



avg_ej = sum_e/time;
plot(1:5000, 10*log10(avg_ej), 'k');
legend('JPPSOV');

% save jppsov0.4.mat
