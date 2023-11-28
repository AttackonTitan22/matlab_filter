% rand函数及单位均匀分布
n=10000;
x=rand(1,n);  %产生(0-1)单位均匀信号，1行，n列
subplot(2,1,1)
plot(x);  %输出信号图
set(gca,'FontSize',20);
title('0-1服从均匀分布的随机序列信号');
subplot(2,1,2)
histogram(x,50)
set(gca,'FontSize',20);
title('0-1服从均匀分布的随机序列直方图');
