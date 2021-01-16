
function [dbmis_GLD,b,C] = GLD(mic_in,noise,m,n_d,f,maxlag,F,C,NFFT,number,k,n_j_i)

%------ Calculate cross-correlation and auto-correlation
crosscorr = xcorr(mic_in(m-n_d+length(f):m-1),noise(m-n_d+length(f):m-1),maxlag);%crosscorrelation of y(n) and white noise i.e. Ryx
lev_noise_input = noise(m-n_d+length(f):m-1);%white noise after softening
noise_cor = xcorr(lev_noise_input,maxlag);%autocoreelation of white noise i.e. Rxx
yn_corr = crosscorr(ceil(length(crosscorr)/2):end);
C=C+1;

%------- Generalized Levinson Durbin Algorithm
p = length(yn_corr);
E = zeros(1,p);
K = zeros(1,p);
a = zeros(1,p);
b = zeros(1,p);
difference = zeros(p,p);
pre_a = zeros(1,p);
pre_b = zeros(1,p);
rxx = xcorr(lev_noise_input',p);
R = rxx(p+2:end);
R_reverse = wrev(R); %Reverse of R
E0 = rxx(p+1);
for M = 1:p
    if M == 1
        K(M) = -R(M)/E0;
        a(M) = K(M);
        b(M) = yn_corr(M)/E0;
        E(M) = (1-(K(M))^2)*E0;
    else
        K(M) = -(R(M) + pre_a(1:M-1) * R_reverse(end-M+2:end))/E(M-1);
        a(M) = K(M);
        b(M) = (yn_corr(M) - pre_b(1:M-1) * R_reverse(end-M+2:end))/E(M-1);
        a(1:M-1) = pre_a(1:M-1) + K(M) * conj(pre_a(M-1:-1:1));
        b(1:M-1) = pre_b(1:M-1) + b(M) * conj(pre_a(M-1:-1:1));
        E(M) = (1-(K(M))^2) * E(M-1);
    end
    pre_a = a;
    pre_b = b;
    
    % Misalignment calculation for each order
    B_FFT=fft(b,NFFT);
    misalignment_GLD=((norm(F(1:NFFT/2)-B_FFT(1:NFFT/2)))^2)/(norm(F(1:NFFT/2))^2);
    dbmis_GLD(number,M,floor(k/n_j_i)+1)=pow2db(misalignment_GLD+10^(-50));
end
end