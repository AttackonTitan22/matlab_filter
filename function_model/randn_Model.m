
% randn生成标准正态分布伪随机数
y=randn(1,10000); 
 subplot(2,1,1);plot(y); 
 set(gca,'FontSize',20);
title('服从高斯分布的随机序列信号'); 
 subplot (2,1,2);histogram(y,50); 
 set(gca,'FontSize',20);
title('服从高斯分布的随机序列直方图');
