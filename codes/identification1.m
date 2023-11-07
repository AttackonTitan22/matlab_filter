%%%%PNSF&SAF&JPPSOV在有色(AR)信号输入非线性系统系统辨识中的应用:Experiment 1
clc;
clear;


deltax=0.2;%插值间距
Slope=1;%初始化斜率
gain=2;%区域边界限制
afinit=0;
table_length=TabFuncLen(deltax,gain,Slope,afinit);%查找表长度
C=1/6*[-1 3 -3 1;3 -6 3 0;-3 0 3 0;1 4 1 0];%B样条基矩阵
% C=1/2*[-1 3 -3 1;2 -5 4 -1;-1 0 1 0;0 2 0 0];%CR样条基矩阵

runs=300;%实验仿真次数
N=5000;%迭代次数
snr=30;%输出信噪比

epsi=1;





L1=14;%NSF线性滤波器长度

% M1=3;%PNSF模块数
% p=1;%外部输入长度
% L2=4;%FIR滤波器长度


M=5;%PNSF模块数
L=3;%模块外部输入数
mw=0.05;
mq=0.5;
mh=0.9;%初始化PNSF参数 ，colored signal

mw_NSF=0.5;
mq_NSF=0.3;%初始化NSF参数,colored signal

p=3;%外部输入信号
L2=1+(p+1)+(p+1)*(p+2)/2;%JPPSOV每个模块经过SOV扩展后的长度
M1=5;%RPSOV模块数
lante1=0.15;%非线性
lante_linear=0.08;%线性

%% PNIIR parameters setting
   p2=3;
   intal_non1=0.008; %0.2;  %27.8483
   intal_linear1=0.01;%0.0008;
   intal_non2=0.008;   %0.23 29.2636； 0.26 27.8790
   intal_linear2=0.01;


if N<30000      % Batch for evaluating MSE
        B=100;
else
        B=4000;
end

em_hpnsf=zeros(N,1);
em_nsf=zeros(N,1);
em_rpsov=zeros(N,1);
em_pniir1=zeros(N,1);
em_pniir2=zeros(N,1);



for run=1:runs
    run
    %%初始化
%     x=zeros(N,1);
    d=zeros(N,1);
   
    
  %% PNSF 
    w_hpnsf=zeros(L+1,M);%每个模块的FIR滤波器的权值不同
    w_hpnsf(1,:)=1.0;
    q=zeros(table_length,M);
    s=zeros(M,1);
    z=zeros(M,N);%每个模块每个时刻的输出
    H_hpnsf=0.001*rand(M,1);    
    y_hpnsf=zeros(N,1);%初始化滤波器输出
    e_hpnsf=zeros(N,1);%初始化误差
 %% SAF  
    w_nsf=zeros(L1,1);
    w_nsf(1)=1;
    q_nsf=zeros(table_length,1);
    s_nsf=zeros(N,1);
    y_nsf=zeros(N,1);
    e_nsf=zeros(N,1);
    
 %% JPPSOV  
    y_rpsov=zeros(N,1);
    e_rpsov=zeros(N,1);w_rpsov=0.001*randn(M1,1);h_rpsov=0.001*randn(L2,1);
    
 %% PNIIR
    r1=zeros(1,M);                            %反馈输入信号(1)：前一个模块的输出
    y_hd1=zeros(M,N);                         %HCNFIR(1)每个模块的非线性输出信号
    h_hd1=zeros(M,p2+2);                       %非线性部分的权值初始化
    w_hd1=zeros(1,M);                         %线性部分的权值初始化
    input_module_hd1=zeros(1,p2+2);            %模块输入信号的初始化
    y_w_hd1=zeros(N,1);
    e_hd1=zeros(N,1);
    x_hd1=zeros(M,p2+2);x_pd1=zeros(M,p2+2);
    yfinal_hd1=zeros(N,1);
    
    r2=zeros(1,M);                            %反馈输入信号(1)：前一个模块的输出
    y_hd2=zeros(M,N);                         %HCNFIR(2)每个模块的非线性输出信号
    h_hd2=zeros(M,p2+2);                       %对于非线性部分权值初始化为0，（L+1）个输入
    w_hd2=zeros(1,M);                         %对于线性部分——各模块权值向量相同时，权值初始化为0，对 M 个模块的输出进行调整    
    input_module_hd2=zeros(1,p2+2);            %模块输入信号的初始化
    y_w_hd2=zeros(N,1);
    yfinal_hd2=zeros(N,1);
    e_hd2=zeros(N,1);
    x_hd2=zeros(M,p2+2);x_pd2=zeros(M,p2+2);
    %%
     %初始化q矩阵
    LutSlope=(table_length-1)/2.0 ; % New slope2
    XX=-LutSlope*deltax;
   for i=1:table_length % Table_Length
    q(i,:)=FUNC(XX,gain,Slope,afinit);
    q_nsf(i)=FUNC(XX,gain,Slope,afinit);
%     qq(i,:)=FUNC(XX,gain,Slope,afinit);
    XX=XX+deltax;
   end
   %%
%    期望响应
   x=randn(N,1);
%      xa=randn(N,1);
%      x=filter(sqrt(1-0.95^2),[1,-0.95],xa);
%      st=filter([0.0154,0.0462,0.0462,0.0154],[1,-1.99,1.572,-0.4583],x);
     
     x=filter(1,[1,-1.79,1.85,-1.27,0.41],x); 
%      x=rand(N,1);
     x=0+(x-min(x))*(0.5-0)/(max(x)-min(x));%将输入向量归一化到[0,0.1]区间
%     x=-0.1+(x-min(x))*(0.1-(-0.1))/(max(x)-min(x));%将输入向量归一化到[-0.1,0.1]区间
    for i=2:N
        d(i)=(d(i-1)/(1+d(i-1).^2))+x(i).^3;       
    end
%     d=sin(st);
%     d=awgn(d,snr,'measured');
%       d=d+dn;
     signy_power=mean(d.^2);
     bsigma=sqrt(signy_power*10^(-snr/10));
     bn=bsigma*randn(N,1);
     d=d+bn;%期望响应
    %%
    for n=L1+M+2:N  %ii代表第ii时刻
           %% PNSF
            u=zeros(M,1);
            UU=zeros(4,M);%由u组成的矩阵，包括4行M列，每列代表一个模块的向量
            UU1=zeros(4,M);
            indexi=zeros(M,1);
            r=zeros(M,1);
            zii=[];
            X=zeros(L+1,M);
            XX1=zeros(L+1,1);
            z1=zeros(M,1);
            hz1x=zeros(L+1,1);
            r_PSOV=zeros(M1,1);y1=zeros(M1,1);temp_w_sum=zeros(L2,1);
            
        for j=1:1:M   %j代表第j个模块
         
           if j==M
               r(j)=x(n-j-1);   %PNSF的另一个输入
           else
               r(j)=z(j+1,n-1); %一步延时
           end
             %%%%%%%%PNSF   
             X(:,j)=[x(n-j:-1:n-j-L+1);r(j)];
             XX1=XX1+X(:,j);
             s(j)=w_hpnsf(:,j)'*X(:,j);%每个模块的线性输出部分
             u(j)=s(j)/deltax-floor(s(j)/deltax);
             indexi(j)=floor(s(j)/deltax)+(table_length-1)/2;
             if indexi(j)<1       
                 indexi(j)=1;     %索引范围必须从1开始
             end
             if indexi(j)>(table_length-3)
                indexi(j)=table_length-3; %索引不能超过Q-3
             end 
            UU(:,j)=[u(j)^3 u(j)^2 u(j) 1]';%第j个模块的u组成的向量
            UU1(:,j)=[3*u(j)^2 2*u(j) 1 0]';%对UU求偏导得到的向量
            z(j,n)=UU(:,j)'*C*q(indexi(j):(indexi(j)+3),j);
            z1(j)=UU1(:,j)'*C*q(indexi(j):indexi(j)+3,j);
%             hz1x=hz1x+H_hpnsf(j)*z1(j)*X(:,j);
%             zii=[zii z(j,n)];
        end
%          zii=zii';
         y_hpnsf(n)=H_hpnsf'*z(:,n)+sum(s);%ii时刻滤波器的输出 
         e_hpnsf(n)=d(n)-y_hpnsf(n);
    %PNSF权值更新
    for i=1:1:M
       
        q(indexi(i):indexi(i)+3,i)=q(indexi(i):indexi(i)+3,i)+mq*e_hpnsf(n)*H_hpnsf(i)*C'*UU(:,i);
        w_hpnsf(:,i)=w_hpnsf(:,i)+mw*e_hpnsf(n)*( H_hpnsf(i)*z1(i)* X(:,i)*(1/deltax)+ X(:,i) );
    end
     H_hpnsf=H_hpnsf+mh*e_hpnsf(n)*z(:,n);
%      w_hpnsf=w_hpnsf+mw*e_hpnsf(n)*(XX1+hz1x);
     
      %% SAF
      s_nsf(n)=w_nsf'*x(n:-1:n-L1+1);
      u1=s_nsf(n)/deltax-floor(s_nsf(n)/deltax);
      tempi_nsf=floor(s_nsf(n)/deltax)+(table_length-1)/2;
            if tempi_nsf<1       
              tempi_nsf=1;     %索引范围必须从1开始
            end
            if tempi_nsf>(table_length-3)
              tempi_nsf=table_length-3; %索引不能超过Q-3
            end      
      y_nsf(n)=[u1^3 u1^2 u1 1]*C*(q_nsf(tempi_nsf:1:tempi_nsf+3));%ii时刻滤波器的输出
      e_nsf(n)=d(n-1)-y_nsf(n);
    %权值更新
    w_nsf=w_nsf+mw_NSF*e_nsf(n)*[3*u1^2 2*u1 1 0]*C*(q_nsf(tempi_nsf:1:tempi_nsf+3))*x(n:-1:n-L1+1);%权值更新
    q_nsf(tempi_nsf:1:tempi_nsf+3)=q_nsf(tempi_nsf:1:tempi_nsf+3)+mq_NSF*e_nsf(n)*C'*[u1^3 u1^2 u1 1]';%控制点更新
   %% JPPSOV 
    for j=1:1:M1   %j代表第j个模块
         
           if j==M1
%                r(j)=z(j,ii-1);   %PNSF的另一个输入
               r_PSOV(j)=x(n-j-1);    %PSOV的另一个输入
%                r_PSOV(j)=y(j);
           else
%                r(j)=z(j+1,ii);
               r_PSOV(j)=y1(j+1);
%              r_PSOV(j)=y(j+1);
           end 
            external_s=[r_PSOV(j);x(n-j:-1:n-(j+p-1))]';
              %%%%%%%%%%% Volterra expansion  
                Xpre=external_s;
                Lpsov=p+1;
                Xpost1=[Xpre];
                  for k=1:Lpsov
                     for kk=k:Lpsov
                         Xpost1=[Xpost1 Xpre(k)*Xpre(kk)];
                     end
                  end  
               new_x=Xpost1;
               XN=[1 new_x].';%经过SOV扩展后的第j个模块的输入向量
               y1(j)=h_rpsov'*XN; %第j个模块的输出
              temp_w_sum=temp_w_sum+w_rpsov(j)*XN;
    end
    
      %%%%%%PSOV输出 
         y_rpsov(n)=w_rpsov'*y1;
         e_rpsov(n)= d(n)-y_rpsov(n);
         w_rpsov=w_rpsov+lante_linear*e_rpsov(n)*y1;  %线性部分更新  

         %   h=h+2*lante1*weight_sum; 
         h_rpsov=h_rpsov+lante1*e_rpsov(n)*temp_w_sum;  %非线性部分更新
         
      %% PNIIR
      %%%%%%%%%%%% 非线性部分，每个模块的计算 %%%%%%%%%%%%%%%%%%%%%    
       g1=zeros(M,1);g2=zeros(M,1);
           for j=M:-1:1
                 g1(j)=y_hd1(1,n-j);
                 g2(j)=y_hd2(j,n-1);
              if j==M
                 r1(j)=y_hd1(1,n-j-1);
                 r2(j)=y_hd2(j,n-2);
              else
                 r1(j)=y_hd1(j+1,n);
                 r2(j)=y_hd2(j+1,n);
              end
      %%%%%%%%%%% PNIIR（1）-外部信号为g1(j) %%%%%%%%%%%%%%%%%%%%
             input_module_hd1=[x(n-j:-1:n-j-p2+1); r1(j); g1(j)]'; 
             x_hd1(j,:)= input_module_hd1;                
             z_hd1(j,n)=x_hd1(j,:)*h_hd1(j,:)';
             y_hd1(j,n)=1./(1+exp(-z_hd1(j,n)));
       
             dif_hd1=x_hd1(j,:)*exp(-z_hd1(j,n))/((1+exp(-z_hd1(j,n))).^2);      %对y_AN1(j,i)求导数

      %%%%%%%%%%% PNIIR（2）-外部信号为p1(j) %%%%%%%%%%%%%%%%   
             input_module_hd2=[x(n-j:-1:n-j-p2+1); r2(j); g2(j)]'; 
             x_hd2(j,:)= input_module_hd2;                
             z_hd2(j,n)=x_hd2(j,:)*h_hd2(j,:)';                                  %第 j 个模块的输出信号(每个模块的输入信号和权值相乘)
             y_hd2(j,n)=1./(1+exp(-z_hd2(j,n)));                                 %第 j 个模块的非线性输出信号
          
             dif_hd2=x_hd2(j,:)*exp(-z_hd2(j,n))/((1+exp(-z_hd2(j,n))).^2);      %对y_AN1(j,i)求导数

          end


%%%%%%%%%%%%%%%%%%%%%%% 线性部分，每个模块的计算 误差计算 ：对于 i 时刻的线性输出的计算 %%%%%%%%%%%%%%%%
      %%%%%%% PNIIR(1)-NLMS-H(n)/W(n) %%%%%%%%%%%%%%%% 
            y_w_hd1(n)=w_hd1*y_hd1(:,n);
            yfinal_hd1(n)=y_w_hd1(n)+sum(z_hd1(:,n));    %线性部分和非线性部分的输出总和
            e_hd1(n)=d(n-1)-yfinal_hd1(n);
   
%     for hi=1:M
%          h_hd1(hi,:)=h_hd1(hi,:)+intal_non1*(1+dif_hd1*w_hd1(hi))*x_hd1(hi,:)*e_hd1(i);  %NLMS算法
%     end  

            for hi=1:M                                                
                x_pd1(hi,:)=x_hd1(hi,:)+dif_hd1*w_hd1(hi);
                h_hd1(hi,:)=h_hd1(hi,:)+intal_non1*e_hd1(n)*x_pd1(hi,:);  %非线性部分——不同模块的权值更新公式（NGD算法）
            end  

            w_hd1=w_hd1+intal_linear1*e_hd1(n)* y_hd1(:,n)';    
    
           

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PNIIR(2)-NLMS-H(n)/W(n) %%%%%%%%%%%%%
            y_w_hd2(n)=w_hd2*y_hd2(:,n);
            
            yfinal_hd2(n)=y_w_hd2(n)+sum(z_hd2(:,n));    %线性部分和非线性部分的输出总和  
            e_hd2(n)=d(n)- yfinal_hd2(n);
   
%     for hi=1:M
%          h_hd2(hi,:)=h_hd2(hi,:)+intal_non1*(1+dif_hd2*w_hd2(hi))*x_hd2(hi,:)*e_hd2(i);  %NLMS算法
%     end  

           for hi=1:M                                               
               x_pd2(hi,:)=x_hd2(hi,:)+dif_hd2*w_hd2(hi);
               h_hd2(hi,:)=h_hd2(hi,:)+intal_non2*e_hd2(n)*x_pd2(hi,:);  %非线性部分——不同模块的权值更新公式（NGD算法）
           end  

           w_hd2=w_hd2+intal_linear2*e_hd2(n)* y_hd2(:,n)';   
           
    end

    
               
   %%
    em_hpnsf=em_hpnsf+e_hpnsf.^2;
    em_nsf=em_nsf+e_nsf.^2;
    em_rpsov=em_rpsov+e_rpsov.^2;
    em_pniir1=em_pniir1+e_hd1.^2;
    em_pniir2=em_pniir2+e_hd2.^2;
    
end
em_hpnsf=em_hpnsf/runs;
em_nsf=em_nsf/runs;
em_rpsov=em_rpsov/runs;
em_pniir1=em_pniir1/runs;
em_pniir2=em_pniir2/runs;

mse_hpnsf=mean(em_hpnsf(end-B-M-1:end-M-1));% Average MSE
mse_nsf=mean(em_nsf(end-B-M-1:end-M-1));% Average MSE
% mse_PNSF=mean(average_e_PNSF(end-B-M-1:end-M-1));% Average MSE
mse_rpsov=mean(em_rpsov(end-B-M-1:end-M-1));% Average MSE
mse_pniir1=mean(em_pniir1(end-B-M-1:end-M-1));% Average MSE
mse_pniir2=mean(em_pniir2(end-B-M-1:end-M-1));% Average MSE

fprintf('PNSF with x range=[0,0.1] Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse_hpnsf,10*log10(mse_hpnsf));
fprintf('SAF with x range=[0,0.1] Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse_nsf,10*log10(mse_nsf));
fprintf('JPPSOV with x range=[0,0.1] Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse_rpsov,10*log10(mse_rpsov));
fprintf('PNIIR(1) with x range=[0,0.1] Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse_pniir1,10*log10(mse_pniir1));
fprintf('PNIIR(2) with x range=[0,0.1] Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse_pniir2,10*log10(mse_pniir2));
% fprintf('PSOV with x range=[0,0.1] Steady-state MSE = %5.7f, equal to %5.3f dB\n',mse_PNSF,10*log10(mse_PNSF));

tmse_hpnsf=10*log10(em_hpnsf);
tmse_nsf=10*log10(em_nsf);
tmse_rpsov=10*log10(em_rpsov);
tmse_pniir1=10*log10(em_pniir1);
tmse_pniir2=10*log10(em_pniir2);

% tmse_PNSF=10*log10(average_e_PNSF);

% plot(1:1:N,10*log10(em_HPNSF),'-r');
% hold on
% plot(1:1:N,10*log10(em_NSF),'-b');
% hold on
% plot(1:1:N,10*log10(average_e_PNSF),'-k');
% legend('HPNSF','NSF','PNSF');
% length=100;
% for i=101:1:N-50
%     tmse_HPNSF(i)=mean(tmse_HPNSF(i-50):1:tmse_HPNSF(i+50));
%     tmse_NSF(i)=mean(tmse_NSF(i-50):1:tmse_NSF(i+50));
%     tmse_PNSF(i)=mean(tmse_PNSF(i-50):1:tmse_PNSF(i+50));
% end
% plot(1:1:N,tmse_HPNSF,'-r');
% hold on
% plot(1:1:N,tmse_NSF,'-b');
% hold on
% plot(1:1:N,tmse_PNSF,'-k');
% % legend('HPNSF','NSF','PNSF');
figure
plot(1:1:N,tmse_hpnsf,1:1:N,tmse_nsf,1:1:N,tmse_rpsov,1:1:N,tmse_pniir1,1:1:N,tmse_pniir2);
 legend('PNSF','SAF','JPPSOV','PNIIR(1)','PNIIR(2)');