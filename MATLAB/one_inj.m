clc
clear all;

%reading input audio file
[ input , Fs ]=audioread( 'alan_1.wav' );
input=input/norm(input);


l = length(input);
n = 20*1e-3*Fs;         % 20ms frames, frame length
n_F = floor(l/n);       % Number of frames


%loading actual feedback path model
load starkeymodelBTE.mat
f=y;
NFFT=1024;
F=fft(f,NFFT);

%specifying the noise durstion and noise injection interval
noise_duration = 40;
noise_injection_interval = 4000;
n_d = Fs*noise_duration/1000;
n_j_i = Fs*noise_injection_interval/1000;
ratio=n_d/n_j_i;

%reading white noise to be injected
load 'input_noise2';

% Make the noise zero mean
for j=0:round(length(r)/n_j_i)-1
    r(j*n_j_i+1:j*n_j_i+n_d)= r(j*n_j_i+1:j*n_j_i+n_d)-mean( r(j*n_j_i+1:j*n_j_i+n_d));
    %     j*n_j_i+1,
    %     j*n_j_i+n_d,
end

%gain values
gain=15; %Forward gain
gain_delay=80; % Forward delay

%constants
m=1;
C=1;
av_num = 1;
soft_idx = 1;

fil = 100; % Filter length
maxlag=fil-1;
dbmis_GLD = zeros(1,fil);

noise_cor = zeros(1,2*maxlag+1);
yn_corr = zeros(1,maxlag+1);
rxx = zeros(2*maxlag+3,1);

% Vectors
u=zeros(1,n);%-----------------------Loudespeaker signal
v=zeros(1,n);%-----------------------I/p for feedback canceller
y=zeros(1,n);%-----------------------Microphone input
x_hat=zeros(1,n);%-------------------Error (or the forward path signal)
u_f=zeros(1,n);%---------------------Feedback signal
z=zeros(1,n);%-----------------------Feedback Canceller O/p
noise=zeros(1,ceil(ratio*length(input)));
mic_in=zeros(1,ceil(ratio*length(input)));

%------- Buffers
U_f=zeros(length(f),1);
V_f=zeros(fil,1);

%------ Soft switching between noise and speech
soft_prt1 = exp(0.05:.05:5)/exp(5);
soft_prt2 = ones(1,n_d-2*length(soft_prt1));
soft_prt3 = wrev(soft_prt1);
soft = cat(2,soft_prt1,soft_prt2,soft_prt3);

number=1;

for i = 1
    F_start(i) = (i-1)*n+1;
    F_end(i) = i*n;
    Frame =input(F_start(i):F_end(i));
    i
    for k=1:length(Frame)% k =sample no. of first frame
        
        if (k==640)
            k=k;
        end
        
        if (mod(k,(C-1)*n_j_i)<=n_d) && (k<=(av_num-1)*n_j_i+n_d)
            u(k)=r(k);%----------------------------------loudespeaker signal = white noise
            % Convolution
            U_f=[u(k); U_f(1:length(f)-1)];
            u_f(k)=f*U_f;%-------------------------------Feedback signal
            y(k)=Frame(k)+u_f(k);%-----------------------Microphone input
            x_hat(k)=y(k);
            
            if mod(k,n_j_i)>length(f)
                noise(m)=r(k)*soft(soft_idx);
                mic_in(m)=y(k);
                m=m+1;
                soft_idx = soft_idx + 1;
            else
                soft_idx = 1;
            end
            
            if mod(k,(C-1)*n_j_i)==n_d
                [dbmis_GLD,b] = GLD(mic_in,noise,m,n_d,f,maxlag,F,C,NFFT,number,k,n_j_i);
            end
            
        else
            v(k)=gain*x_hat(k-gain_delay);
            u(k)=v(k);
            U_f=[u(k); U_f(1:length(f)-1)];
            u_f(k)=f*U_f;
            y(k)=Frame(k)+u_f(k);
            
            V_f=[v(k); V_f(1:length(b)-1)];
            z(k)=b*V_f;
            x_hat(k)=y(k)-z(k);
        end
    end
end

prev_U_f = U_f;
prev_x_hat = x_hat;
final_x_hat = x_hat;
prev_V_f = V_f;
final_noise = noise(1:n-length(f));
final_mic_in = mic_in(1:n-length(f));
hh=m;
m=1;

for i = 2:n_F
    F_start(i) = (i-1)*n+1;
    F_end(i) = i*n;
    Frame =input(F_start(i):F_end(i));
    i
    
    if (i==n_F)
        i=i;
    end
    
    for k=1:length(Frame) % k =sample no. of first frame
        
        if (k==81)
            k=k;
        end
        
        if (mod(((i-1)*n+k),(C-1)*n_j_i)<=n_d) && (((i-1)*n+k)<=(av_num-1)*n_j_i+n_d)
            u(k)=r(((i-1)*n+k));%----------------------------------loudespeaker signal = white noise
            % Convolution
            U_f=[u(k); prev_U_f(1:length(f)-1)];
            u_f(k)=f*U_f;%-------------------------------Feedback signal
            y(k)=Frame(k)+u_f(k);%-----------------------Microphone input
            x_hat(k)=y(k);
            
            
            if mod(((i-1)*n+k),n_j_i)>length(f)
                noise(m)=r(((i-1)*n+k))*soft(soft_idx);
                mic_in(m)=y(k);
                m=m+1;
                hh=hh+1;
                soft_idx = soft_idx + 1;
            else
                soft_idx = 1;
            end
            
            
            
            if mod(((i-1)*n+k),(C-1)*n_j_i)==n_d
                final_noise = [final_noise,noise];
                final_mic_in = [final_mic_in,mic_in];
                [dbmis_GLD,b] = GLD(final_mic_in,final_noise,hh,n_d,f,maxlag,F,C,NFFT,number,k,n_j_i);
            end
        else
            
            if (k-gain_delay<=0)
                v(k)=gain*prev_x_hat(n+k-gain_delay);
            else
                v(k)=gain*prev_x_hat(k-gain_delay);
            end
            u(k)=v(k);
            U_f=[u(k); prev_U_f(1:length(f)-1)];
            u_f(k)=f*U_f;
            y(k)=Frame(k)+u_f(k);
            
            V_f=[v(k); prev_V_f(1:length(b)-1)];
            z(k)=b*V_f;
            x_hat(k)=y(k)-z(k);
        end
        
        prev_U_f = U_f;
        prev_x_hat = x_hat;
        prev_V_f = V_f;
        
    end
    
    
    final_x_hat = [final_x_hat,x_hat];
    m=1;
    
end



figure( 1 );
avg_dbmis_GLD = dbmis_GLD;
plot( avg_dbmis_GLD );
xlabel('Filter Length (M)');
ylabel('Misalignment (dB)');
grid on;

py_value = pesq(input,final_x_hat);

% %signal with feedback
% load ('e_feedback.mat');
% figure( 2 );
% subplot(2,1,2);spectrogram(x_hat/norm(x_hat),512,500,512,Fs,'yaxis','reassigned');title('after feedback cancellation');
% subplot(2,1,1);spectrogram(tt/norm(tt),512,500,512,Fs,'yaxis','reassigned');title('with feedback');