function e =PNIIR1_Model(L,M,N,A,y,x,d,lrh,lrw,h,w)

        % L:外部输入个数
        % M:模块数
        % N:信号数
        % A:非线性部分权重更新算子
        % y:每个模块的输出信号
        % x:输入信号
        % d:期望输出
        % lrh:非线性部分的初始权重
        % lrw:线性部分的初始权重
        % h:非线性部分的权重矩阵
        % w:线性部分的权重矩阵

        for i=L+M:1:N        %输入信号x的长度
            Z=0;
            for j=M:-1:1    %模块的长度
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
                 %非线性部分权重更新
                h(j,:)=h(j,:)+lrh*e(i)* A(j,:)/((norm(A(j,:)).^2));  
            end      
                 %线性部分权重更新
            %w=w+lrw*e(i)*y(i,:)/(norm(y(i,:)).^2);
            w=w+lrw*e(i)*y(i,:)/(norm(y(i,:)).^2);
               
        end
end
