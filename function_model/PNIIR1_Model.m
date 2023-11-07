function e =PNIIR1_Model(L,M,N,A,y,x,d,lrh,lrw,h,w)

        % L:�ⲿ�������
        % M:ģ����
        % N:�ź���
        % A:�����Բ���Ȩ�ظ�������
        % y:ÿ��ģ�������ź�
        % x:�����ź�
        % d:�������
        % lrh:�����Բ��ֵĳ�ʼȨ��
        % lrw:���Բ��ֵĳ�ʼȨ��
        % h:�����Բ��ֵ�Ȩ�ؾ���
        % w:���Բ��ֵ�Ȩ�ؾ���

        for i=L+M:1:N        %�����ź�x�ĳ���
            Z=0;
            for j=M:-1:1    %ģ��ĳ���
                if j==M
                    r(j)=y(i-j-1,1);
                    g(j)=y(i-j,1);
                else
                    r(j)=y(i,j+1);
                    g(j)=y(i-j,1);
                end
                X=[x(i-j:-1:i-(j+L-1)) r(j) g(j)]; 
                z(i,j)=h(j,:)*X.';
                y(i,j)=1/(1+exp(-z(i,j)));
%                 syms xn;
%                 f(xn)=1/(1+exp(-xn));
%                 f(xn)=diff(f(xn),xn);
%                 A(j,:)=X+w(j)*f(h(j,:)*X.')*X;
                A(j,:)=X+w(j)*(exp(-h(j,:)*X.')/(exp(-h(j,:)*X.') + 1)^2)*X;
                Z=Z+z(i,j);
            end
            y_final(i)=w*y(i,:).'+Z;
            e(i)=d(i-1)-y_final(i);
            for j=1:1:M
                 %�����Բ���Ȩ�ظ���
                h(j,:)=h(j,:)+lrh*e(i)* A(j,:)/((norm(A(j,:)).^2));  
            end      
                 %���Բ���Ȩ�ظ���
            %w=w+lrw*e(i)*y(i,:)/(norm(y(i,:)).^2);
            w=w+lrw*e(i)*y(i,:)/(norm(y(i,:)).^2);
               
        end
end
