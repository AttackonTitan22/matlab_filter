function d=nonlinear1(x,N)
        d=zeros(1,N);
        for n=2:N
               d(n)=d(n-1)/(1+d(n-1).^2)+x(n).^3;
        end
end