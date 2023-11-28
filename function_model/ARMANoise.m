function x=ARMANoise(N)
        x=zeros(1,N);
        u=randn(1,N);
        x(1)=-0.1*u(1);
        x(2)=0.04*x(1)-0.1*u(2)-0.01*u(1);
        x(3)=0.04*x(2)-0.034*x(1)-0.1*u(3)-0.01*u(2)-0.137*u(1);
        x(4)=0.04*x(3)-0.034*x(2)+0.0396*x(1)-0.1*u(4)-0.01*u(3)-0.137*u(2)+0.0353*u(1);
        x(5)=0.04*x(4)-0.034*x(3)+0.0396*x(2)-0.07565*x(1)-0.1*u(5)-0.01*u(4)-0.137*u(3)+0.0353*u(2)+0.06984*u(1);
        for n=6:N
               x(n)=0.04*x(n-1)-0.034*x(n-2)+0.0396*x(n-3)-0.07565*x(n-4)-0.1*u(n)-0.01*u(n-1)-0.137*u(n-2)+0.0353*u(n-3)+0.06984*u(n-4);
        end
end