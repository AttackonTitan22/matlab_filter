function x=ARMANoise2(N)
        x=zeros(N);
        u=randn(N);
        x(1)=u(1);
        x(2)=1.79*x(1)+u(2);
        x(3)=1.79*x(2)-1.85*x(1)+u(3);
        x(4)=1.79*x(3)-1.85*x(2)+1.27*x(1)+u(4);
        for n=5:N
               x(n)=1.79*x(n-1)-1.85*x(n-2)+1.27*x(n-3)-0.41*x(n-4)+u(n);
        end
        
end