%运用于PNFIR采用EFAST方法进行敏感性分析
function d=non2(x,N)
        d=zeros(1,N);
        for n=3:N
               d(n)=d(n-1)/(1+d(n-1).^2)+x(n).^2+x(n-2).^2;
        end
end