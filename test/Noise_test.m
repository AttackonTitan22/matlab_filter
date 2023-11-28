
clc;
clear;
close all;

N=5000;
x=0.6*rand(1,N);
subplot(2,3,1)
plot(x);  %输出信号图
set(gca,'FontSize',8);
title('0-0.6服从均匀分布的随机序列信号');
subplot(2,3,4)
hist(x,50)
set(gca,'FontSize',8);
title('0-0.6服从均匀分布的随机序列直方图');
mean_x = mean(x) %验证均值为0.5  0-1
power_x = var(x) %验证功率为1/12  0-1


xn=randn(1,N);xn=0.6*(xn-min(xn))/(max(xn)-min(xn));
subplot(2,3,2)
plot(xn);  %输出信号图
set(gca,'FontSize',8);
title('0-0.6服从均值为0方差为1的高斯序列信号');
subplot(2,3,5)
hist(xn,50)
set(gca,'FontSize',8);
title('0-0.6服从均值为0方差为1的高斯序列直方图');
power_y = var(xn)     %验证功率
mean_y = mean(xn)     %验证均值


xA=ARMANoise(N);  xA=0.6*(xA-min(xA))/(max(xA)-min(xA));  xA=xA.';
subplot(2,3,3)
plot(xA);  %输出信号图
set(gca,'FontSize',8);
title('0-0.6服从指定有色噪声序列信号');
subplot(2,3,6)
hist(xA,50)
set(gca,'FontSize',8);
title('0-0.6服从指定有色噪声序列直方图');


