

%  IEEE  Trans. Signal Processing, Vol.57,No.1,January 2009.
%  Main program
%  Title: A novel adaptive nonlinear filter-based pipelined feedforward
%  second-order Volterra architecture.
%  Simple Pipelined seond order volterra  joint linear filter by NLMS;
%  Including two sections: linear and nonlinear;
%  Written by Haiquan Zhao (Ph.d) 
%  2011.1.18,

%h=0.001*randn(1,N1);这个东西真的尤其重要，因为没有乘0.001，有时候就会不收敛

        clc;
        clear all;
        close all;

        N=5000;
        sum_e_linear = 0;
        sum_ev = 0;
        sum_ep = 0;

time=100;        
for t=1:time
    t
        %x=rand(1,N);  %w=0.9 h=0.3 x(kk).^3
        x=0.4*randn(1,N);x=(x-min(x))/(max(x)-min(x));  %w=0.5 h=0.3  x(kk).^2 (h过小收敛慢，h过大会有脉冲且稳态误差会增加)
%          x=color_arma_signal_new(N);  %w=0.3 h=0.016  0.5*x(kk).^3
%         % x_a=0.36; 
%          x=(x-min(x))/(max(x)-min(x));

        d=zeros(1,N);
        
        %%%%%%%%%%%%%%%%%%%%%%% Nonlinear system identification

        for kk=2:N
               d(kk)=d(kk-1)/(1+d(kk-1).^2)+x(kk).^3;
        end
   
        snr=20;
        info_power=(norm(d)^2)/N;        %信号的能量            
        snr_lin=10^(snr/10);                         
        noise_power=info_power/snr_lin;   %噪声的方差
        d=d+sqrt(noise_power)*randn(1,N); 
        
        % 学习率
%         lante1=0.8;  %jppsov  w
%         lante2=0.2;  %jppsov  h
%         lante1=0.1;  %jppsov  w
%         lante2=0.1;  %jppsov  h
        lante1=0.5;  %jppsov  w
        lante2=0.3;  %jppsov  h
        
        lante3=0.01;  %sov
        lante4=0.01;  %psov
        
        %%%%%%%%%%%%%%%%%%%%%%
        %    Parameters Lists:
        %                       train_x:  trainning sample;
        %                  train_number:  trainning number;
        %                 module_number:  module's number;
        %                  input_number:  input signal length 's number;
        %                      desire_d:  desire of output signal through the unknown system; 
        %                      lt:    parameter list;
        %   Output results variables:
        %                    e_d:  Error of output; 
        %                y_final:  final output  through the processing to the data;
        %                       

        train_x=x;
        train_number=N;
        desire_d=d;
        module_number=5;
        lt1=lante1;
        lt2=lante2;


        % Number of modules
        % M=10;
        M=module_number;
        % Number of each module's input signed
        p=1;
        % Length of trainning sampling.
        %n=5000;
        n=train_number;

        % Output of each module
        y=zeros(1,M);
        y1=zeros(1,M);
        % Length of each module allover input
        N1=1+(p+1)+(p+1)*(p+2)/2;%jppsov

        % Weight Vector matrix of nonlinear filter
        h=0.001*randn(1,N1); %jppsov
        
        hp=0.001*randn(1,N1);%psov
        % external signal of each module
        % q=5;
        % Weight vector of linear section;
        
        w=0.001*randn(1,M);

        external_s=zeros(1,p+1);

%         s=zeros(1,N);

        r=zeros(1,M);


        x=train_x;
        d=desire_d;

        %%%%%%%%%%%%%%%%%%%%%%%


%         lante_linear=lt1;
%         lante1=lt2;
        
        
        
        m1=10;% input_length;sov
        NN=1+m1+m1*(m1+1)/2;
        hh=0.001*randn(1,NN);


        %%%%sov的非线性
        for i=m1:1:N %输入信号x的长度

            %--------Nonlinear subset-------------

            input_xv=[x(i:-1:i-m1+1)];

            
                    %%%%%%%%%%% Volterra expansion  
                    Xpre=input_xv;%[ Xpre(1) Xpre(2) ... Xpre(L)]
                    L=m1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    %     
                    %    对矢量Xpre进行非线性扩展
                    %    最后生成的矢量Xpost为[ 1 Xpre(1) Xpre(2) ... Xpre(L) Xpre(1)*Xpre(1) Xpre(1)*Xpre(2) ... Xpre(L)*Xpre(L)]
                    %
                    % 		
                    % 		Xpre and Xpost are vectors
                    % 		L is the length
                    %

                    %Xpost=[1 Xpre];

                    Xpost=[Xpre];
                    
                    for ii=1:L
                       for jj=ii:L
                          Xpost=[Xpost Xpre(ii)*Xpre(jj)];
                       end
                    end  
                    new_xv=Xpost;
            
         
            sv=[1 new_xv].';  %最终得到的volterra扩展的输入值[ 1 Xpre(1) Xpre(2) ... Xpre(L) Xpre(1)*Xpre(1) Xpre(1)*Xpre(2) ... Xpre(L)*Xpre(L)]
            temp_w_sum=zeros(1,N1);


            for j=M:-1:1  %模块的长度

                %---------Desired signal of j th module-------
                %dd(j)=d(i-j);

                %-------- Input of i th module----------------
                %external_s=[x(i-j:-1:i-(j+p-1))];

                if j==M
                    r(j)=y(j);
                else
                    r(j)=y(j+1);
                end

                external_s=[x(i-j:-1:i-(j+p-1)) r(j)];
                    %%%%%%%%%%% Volterra expansion  
                    Xpre=external_s;
                    L=p+1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    %     
                    %    对矢量Xpre进行非线性扩展
                    %    最后生成的矢量Xpost为[ 1 Xpre(1) Xpre(2) ... Xpre(L) Xpre(1)*Xpre(1) Xpre(1)*Xpre(2) ... Xpre(L)*Xpre(L)]
                    %
                    % 		
                    % 		Xpre and Xpost are vectors
                    % 		L is the length
                    %
                    %Xpost=[1 Xpre];
                    Xpost=[Xpre];
                    for ii=1:L
                       for jj=ii:L
                          Xpost=[Xpost Xpre(ii)*Xpre(jj)];
                       end
                    end  
                    new_x=Xpost;
                    %%%%%%%%%%%%%%%%%%%%%%%
                s=[1 new_x].'; % [ 1 Xpre(1) Xpre(2) ... Xpre(L) Xpre(1)*Xpre(1) Xpre(1)*Xpre(2) ... Xpre(L)*Xpre(L)]
                y(j)=h*s;   %式子3-35
                
                yp(j)=hp*s; 
                
                temp_w_sum=temp_w_sum+w(j)*s.';    

            end
            
            
            y_final(i)=w*y.';
            d_linear(i)=d(i-1);
            e_linear(i)=d_linear(i)-y_final(i);  %式子3-40
            w=w+lante1*e_linear(i)*y/(1+norm(y).^2);    

         %   h=h+2*lante1*weight_sum; 
            h=h+lante2*e_linear(i)* temp_w_sum/((1+norm(temp_w_sum).^2));  
            
            yv(i)=hh*sv;%是sov  sv就是输入扩展信号
            ev(i)=d(i)-yv(i);
            hh=hh+lante3*ev(i)*sv.'/(1+norm(sv).^2);
            wv(i)=norm(hh);
            
            temp_hp_input=s;  %是PSOV
            ep(i)=d(i-1)-yp(1);
            hp=hp+lante4*ep(i)*temp_hp_input.'/((1+norm(temp_hp_input).^2));  
            
        end
        sum_e_linear = sum_e_linear + e_linear.^2;
        sum_ev = sum_ev + ev.^2;
        sum_ep = sum_ep + ep.^2; 
        
end

avg_e_linear = sum_e_linear/time;
avg_ev = sum_ev/time;
avg_ep = sum_ep/time;

% avg_e_linear(1:m1) = 1e-1;
% edb = 10*log10(avg_e_linear);
% emf = zeros(size(edb));
% lambda = 0.998;
% for n = 2:length(emf)
%     emf(n) = lambda*emf(n-1) + (1-lambda)*edb(n);
% end
% 
% avg_ev(1:m1) = 1e-1;
% edb_1 = 10*log10(avg_ev);
% emf_1 = zeros(size(edb_1));
% lambda = 0.998;
% for n = 2:length(emf_1)
%     emf_1(n) = lambda*emf_1(n-1) + (1-lambda)*edb_1(n);
% end
% 
% avg_ep(1:m1) = 1e-1;
% edb_2 = 10*log10(avg_ep);
% emf_2 = zeros(size(edb_2));
% lambda = 0.998;
% for n = 2:length(emf_2)
%     emf_2(n) = lambda*emf_2(n-1) + (1-lambda)*edb_2(n);
% end
        
%%%%%%%%%%%%%%%%%%%%%%

% figure(1);
% plot(1:5000,10*log10(avg_e_linear.^2));  
% figure(2);
% plot(1:5000,10*log10(avg_ev.^2)); 
% figure(3);
% plot(1:5000,10*log10(avg_ep.^2)); 
%         
% figure(1);
% plot(emf);  
% % figure(2);
% % plot(d);
% figure(3);
% plot(emf_1);
% figure(4);
% plot(emf_2);
% figure(3);
% plot(y_final);
% figure(5);
% plot(10*log10(Mse_sum_result));
% figure(1); 
% %
% plot(1:5000, emf, 'g', 1:5000, emf_1, 'b',1:5000, emf_2, 'k');
% legend('JPSOV','SOV','PSOV');
%
plot( 1:5000, 10*log10(avg_ev), 'm--',1:5000, 10*log10(avg_ep), 'c',1:5000, 10*log10(avg_e_linear), 'k');
legend('SOV','PSOV','JPPSOV');
% plot(1:5000, 10*log10(avg_e_linear), 'g');
% legend('JPPSOV');
% plot( 1:5000, 10*log10(avg_ev), 'b',1:5000, 10*log10(avg_ep), 'r');
% legend('SOV','PSOV');
