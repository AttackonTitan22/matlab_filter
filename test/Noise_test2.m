xn=randn(1,N);xn=0.6*(xn-min(xn))/(max(xn)-min(xn));
subplot(2,2,1)
plot(xn);  %输出信号图
set(gca,'FontSize',8);
title('服从均值为0方差为1的高斯序列信号');
subplot(2,2,3)
hist(xn,50)
set(gca,'FontSize',8);
title('服从均值为0方差为1的高斯序列直方图');
power_y = var(xn)     %验证功率
mean_y = mean(xn)     %验证均值


xA=ARMANoise(N);  xA=0.6*(xA-min(xA))/(max(xA)-min(xA));  xA=xA.';
subplot(2,2,2)
plot(xA);  %输出信号图
set(gca,'FontSize',8);
title('服从指定有色噪声序列信号');
subplot(2,2,4)
hist(xA,50)
set(gca,'FontSize',8);
title('服从指定有色噪声序列直方图');
power_y = var(xA)     %验证功率
mean_y = mean(xA)     %验证均值
