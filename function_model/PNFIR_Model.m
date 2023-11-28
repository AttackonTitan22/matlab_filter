function e=PNFIR_Model(L,M,N,A,y,x,d,lrh,lrw,lra,h,w,TT,Ta)

        % L:外部输入个数
        % M:模块数
        % N:信号数
        % A:非线性部分权重更新算子
        % y:每个模块的输出信号
        % x:输入信号
        % d:期望输出
        % lrh:非线性部分的初始权重
        % lrw:线性部分的初始权重
        % lra:凸组合系数的初始权重
        % h:非线性部分的权重矩阵
        % w:线性部分的权重矩阵
        % TT:凸组合系数lambda
        % Ta:凸组合系数lambda的初始值

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
