        clc;
        clear;
        close all;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                               %
        %         x=0.6*rand(1,N);                      %
        %         learning_rate_w=0.01;  %  w           %
        %         learning_rate_h=0.6;   %  h           %
        %         learning_rate_a=0.6;   %  a           %
        %                                               %
        %          ��ʼh=0.001*randn(M,L+1);            %
        %          ��ʼw=0.001*randn(1,M);              %
        %                                               %
        %              TT=1;                            %
        %              Ta=0;                            %
        %                                               %
        %                                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        x=0.6*rand(1,N);  %w=0.01 h=0.6 a=0.6
        %x=randn(1,N);x=0.4*(x-min(x))/(max(x)-min(x)); %w=0.5 h=0.3  x(kk).^2 (h��С��������h���������������̬��������)
        %x=ARMANoise2(N);  x=0.6*(x-min(x))/(max(x)-min(x));  x=x.'; % w=0.001  h=0.5  a=0.01
       
       

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
        %e=PNFIR_Model(L,M,N,A,y,x,d,lrh,lrw,lra,h,w,TT,Ta);      
        for i=L+M-1:1:N        %�����ź�x�ĳ���
            Z=0;
            for j=M:-1:1    %ģ��ĳ���
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
                 %�����Բ���Ȩ�ظ���
                %h(j,:)=h(j,:)+lrh*e(i)* A(j,:)/(1+(norm(A(j,:)).^2)); 
                h(j,:)=h(j,:)+lrh*e(i)* A(j,:); 
            end      
                 %���Բ���Ȩ�ظ���
            w=w+TT*lrw*e(i)*y(i,:); 
                 %͹��ϲ�������
            Ta=Ta+lra*e(i) *(w*y(i,:).'-Z)*TT*(1-TT); 
        end

        sum_e= sum_e + e.^2;
end

avg_e = sum_e/time;
plot(1:5000, 10*log10(avg_e), 'k');
legend('PNFIR');

%save ('PNFIR\data\PNFIR.mat','avg_e');
save ('PNFIR\data\rand\PNFIR0.01_0.6_0.6.mat','avg_e');
%save ('PNFIR\data\PNFIR0.001_0.5_0.01_ARMA.mat','avg_e');