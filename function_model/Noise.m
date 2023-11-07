function d=Noise(snr,d,N)
% snr: 信噪比
% d: 期望输出
% N: 信号数
        info_power=(norm(d)^2)/N;  %信号的能量                  
        snr_lin=10^(snr/10);                         
        noise_power=info_power/snr_lin;   %噪声的方差
        d=d+sqrt(noise_power)*randn(1,N);
end