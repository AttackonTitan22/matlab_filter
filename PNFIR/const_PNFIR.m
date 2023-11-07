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
        external_x=zeros(1,L+1);
        % ÿ��ģ�������ӳ��ź�
        r=zeros(1,M);
        
        for i=L+M-1:1:N        %�����ź�x�ĳ���
            Z=0;
            A=zeros(M,L+1);
            for j=M:-1:1    %ģ��ĳ���
                if j==M
                    r(j)=y(i-1,1);
                else
                    r(j)=y(i-1,j+1);
                end
                external_x=[x(i-j+1:-1:i-(j+L-2)) r(j)]; 
                z(i,j)=h(j,:)*external_x.';
                y(i,j)=1/(1+exp(-z(i,j)));
                %A(j,:)=A(j,:)+w(j)*external_x*(exp(-h(j,:)*external_x.')/(exp(-h(j,:)*external_x.') + 1)^2);
                A(j,:)=w(j)*external_x*(exp(-h(j,:)*external_x.')/(exp(-h(j,:)*external_x.') + 1)^2);

                %Z=Z+z(i,j);
                %TT=1/(1+exp(-Ta));
            end
            y_final(i)=TT*w*y(i,:).'+(1-TT)*Z;
            e(i)=d(i-1)-y_final(i);
            for j=1:1:M
                %�����Բ���Ȩ�ظ���

                %�̶�TTʱ������JPPSOV����
                h(j,:)=h(j,:)+lrh*e(i)* A(j,:)/(1+(norm(A(j,:)).^2)); 
                %�������õĸ����㷨
                %h(j,:)=h(j,:)+lrh*e(i)* A(j,:); 
            end    
                 %���Բ���Ȩ�ظ���
            w=w+TT*lrw*e(i)*y(i,:)/(1+norm(y(i,:)).^2); 
                 %͹��ϲ�������
            %Ta=Ta+lra*e(i)*(w*y(i,:).'-Z)*TT*(1-TT); 
               
        end
        sum_e= sum_e + e.^2;
end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'k');
legend('PNFIR(͹���ϵ��Ϊ1)');
save ('PNFIR\data\const_PNFIR.mat','avg_e');

