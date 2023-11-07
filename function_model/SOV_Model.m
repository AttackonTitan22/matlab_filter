function e=SOV_Model(x,N,M,d,lrh,y,h)
        % x:输入信号
        % M:模块数
        % N:信号个数
        % dedire_d:期望输出
        % lrh:sov非线性部分的权值

        
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

end
