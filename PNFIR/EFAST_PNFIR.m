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
        x=0.5*rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        %x=randn(1,N);x=0.4*(x-min(x))/(max(x)-min(x)); %w=0.5 h=0.3  x(kk).^2 (h过小收敛慢，h过大会有脉冲且稳态误差会增加)
        
        d=non2(x,N);
        d=Noise(30,d,N);
                  

        module_number=5;
        lrw=learning_rate_w;
        lrh=learning_rate_h;
        lra=learning_rate_a;

        % Number of modules
        M=module_number;
        % Length of external input
        L=1;  % L=2
        % Output of each module
        z=zeros(N,M);
        y=zeros(N,M);
        % Weight Vector matrix of nonlinear filter
        % h=0.001*randn(M,L+1);
        h=0.001*randn(M,L+1);
        % Weight vector of linear section;       
        w=0.001*randn(1,M);
        % 每个模块的输入向量
        X=zeros(1,L+1);
        % 每个模块的输出延迟信号
        r=zeros(1,M);
        A=zeros(M,L+1);

        
        for i=L+M-1:1:N        %输入信号x的长度
            Z=0;
            for j=M:-1:1    %模块的长度
                if j==M
                    r(j)=y(i-1,1);
                else
                    r(j)=y(i-1,j+1);
                end
                
                X=[x(i-j+1:-1:i-(j+L-2)) r(j)]; 
                z(i,j)=h(j,:)*X.';
                y(i,j)=1/(1+exp(-z(i,j)));
%                  syms xn;
%                  f(xn)=1/(1+exp(-xn));
%                  f(xn)=diff(f(xn),xn);
%                  A(j,:)=X+w(j)*f(h(j,:)*X.')*X;
                A(j,:)=(TT*w(j)*(exp(-h(j,:)*X.')/(exp(-h(j,:)*X.') + 1)^2)+(1-TT))*X;
                Z=Z+z(i,j);
            end
            TT=1/(1+exp(-Ta));
            y_final(i)=TT*w*y(i,:).'+(1-TT)*Z;
            e(i)=d(i)-y_final(i);
            for j=1:1:M
                 %非线性部分权重更新
                %h(j,:)=h(j,:)+lrh*e(i)* A(j,:)/(1+(norm(A(j,:)).^2)); 
                h(j,:)=h(j,:)+lrh*e(i)* A(j,:); 
            end      
                 %线性部分权重更新
            w=w+TT*lrw*e(i)*y(i,:); 
                 %凸组合参数更新
            Ta=Ta+lra*e(i) *(w*y(i,:).'-Z)*TT*(1-TT); 
        end
        sum_e= sum_e + e.^2;    

end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'k');
legend('PNFIR');

save ('PNFIR\data\PNFIR_L=1.mat','avg_e');