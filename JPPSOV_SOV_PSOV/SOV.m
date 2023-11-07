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
        x=1*randn(1,N);x=(x-min(x))/(max(x)-min(x));  %w=0.5 h=0.3  x(kk).^2 (h��С��������h���������������̬��������)

        d=nonlinear1(x,N);
        d=Noise(20,d,N);
                    
        module_number=10;
        lrh=learning_rate_h;
        % ģ����
        M=module_number;
        % ÿ��ģ������
        y=zeros(1,N);
        % ÿ��ģ����������
        N1=1+M+(M+1)*M/2;
        % �����Բ��ֵ�Ȩ�ؾ���
        h=0.001*randn(1,N1);
        
        for i=M:1:N        %�����ź�x�ĳ���
      
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
            %�����Բ���Ȩ�ظ���
            h=h+lrh*e(i)* X_post.'/((1+norm(X_post).^2));
                      
        end
         sum_e= sum_e + e.^2;
end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'c');
legend('SOV');


