function e=PSOV_Model(N1,N,M,p,x,y,lrh,h,d)


        % N1:每个模块输入信号的个数
        % N:信号个数
        % M:模块数
        % p:每个模块的外部输入个数
        % x:输入信号
        % y:每个模块的输出信号
        % lrh:非线性部分的初始权重
        % h:非线性部分的权重矩阵
        % r:每个模块的输出延迟信号
        % d:期望输出

        for i=N1:1:N        %输入信号x的长度
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

            end
            e(i)=d(i-1)-y(1);

            %非线性部分权重更新
            h=h+lrh*e(i)* X_post.'/((1+norm(X_post).^2));  
        end