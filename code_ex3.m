function code_ex3
c = 340;
%Beam pattern
B =@(theta,N,d,theta_s,f) (1/N) * ((sin((N/2)*((2*pi*f)/c)*d*(cos(theta) - cos(theta_s))))./(sin(((2*pi*f)/(2*c))*d*(cos(theta)-cos(theta_s)))));

%% 1.4 Part 1 
%Plots of magnitude of beam pattern (changing microphone number, d=8 cm)
figure;
for i=4:4:16
    N = i;
    fplot (@(theta) abs(B(theta, N , 0.08, pi/2, 2000)),[0 pi]);
    xlabel('è(rad)');
    ylabel ('Magnitude (dB)');
    hold on;
    set(gca,'YScale','log');
end
suptitle('Delay-and-Sum Beam Pattern');
legend('N = 4' , 'N = 8', 'N = 12', 'N = 16'); 
hold off;
    
%% 1.4 Part2 
%Plots of magnitude of beam pattern (changing distance, N=8) 
figure;
for i = 8 : 4 : 20
    d = i/100;
    fplot (@(theta) abs(B(theta, 8 , d, pi/2, 2000)),[0 pi]);
    xlabel('è (rad)');
    ylabel ('Magnitude (dB)');
    hold on;
    set(gca,'YScale','log');
end
suptitle('Delay-and-Sum Beam Pattern');
legend('d = 8cm' , 'd = 12cm', 'd = 16cm', 'd = 20cm');
hold off;
figure;

%% 1.4 Part3 
%Plots of magnitude of beam pattern (changing theta)
j=1;
for theta = 0 : 15 : 60 
    if (theta >=45) || (theta == 0) 
        fplot (@(f) abs(B(deg2rad(theta), 8 , 0.08, pi/2, f)),[0 8000]);
        xlabel('f (Hz)');
        ylabel ('Magnitude (dB)');
        hold on;
        set(gca,'YScale','log');
        j=j+1;
    end
end
legend('è = 0^ï' ,'è = 45^0' ,'è = 60^ï'); 
hold off; 
suptitle('Delay-and-Sum Beam Pattern');

%% 1.4 Part4 
%Plots of magnitude of beam pattern (changing theta_s)
for theta_s = 0 : 15 : 90 
    if (theta_s ==45) || (theta_s == 0) || (theta_s == 90) 
        figure;
        semilogr_polar((-pi:0.001:pi) , abs( B((-pi:0.001:pi) , 8 , 0.08 , deg2rad(theta_s), 2000))) ;
        legend(['ès =',num2str(theta_s),' ^o']);
        title('Delay-and-Sum Beam Pattern');
        xlabel('è (^o)');
        ylabel('Magnitude (dB)');
    end
end
hold off; 

%% 2.1.A.1
N=7;
c=340;
dist=0.08;
theta_s=pi/4;
fs=48000;
%Audioread of simulated signals
x(:,1)=audioread('sensor_0.wav');
x(:,2)=audioread('sensor_1.wav');
x(:,3)=audioread('sensor_2.wav');
x(:,4)=audioread('sensor_3.wav');
x(:,5)=audioread('sensor_4.wav');
x(:,6)=audioread('sensor_5.wav');
x(:,7)=audioread('sensor_6.wav');
source=audioread('source.wav');
%FFT of beamformer input
X=fft(x);
X=fftshift(X,1);
d_ks=zeros(length(x),7);
w =2*pi*fs*((-length(x)/2:length(x)/2-1)/length(x));
%Calculation of d(ks)
for i=0:N-1
    d_ks(:,i+1)=exp(-1j*((N-1)/2)*(w/c)*dist*cos(theta_s)).*exp(1j*i*(w/c)*dist*cos(theta_s));
end
d_ks=transpose(d_ks);
d_ks=ctranspose(d_ks);
%Frequency response of beamformer
H=(1/N).*d_ks;
Y=H.*X;
Y=ifftshift(Y,1);
y=ifft(Y);
%Beamformer output
y_final=sum(y,2);
y_final=real(y_final);
audiowrite('sim_ds.wav',y_final,48000);

%% 2.A.2
t=(0:length(source)-1)*(1/fs);
%Plot of clean signal
figure;
plot(t,source);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Voice signal without noise");
%Plot of noisy signal
figure;
plot(t,x(:,4));
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Voice signal of microphone 3 with noise");
%Plot of beamformer output
figure;
plot(t,y_final);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Output of delay-and-sum beamformer");
%Spectrogram of clear signal
[~,F1,T1,P1]=spectrogram(source,0.01*fs,0.005*fs,0.01*fs,fs);
figure;
surf(T1.*fs,F1,10*log10(P1),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of voice signal without noise');
%Spectrogram of noisy signal
[~,F2,T2,P2]=spectrogram(x(:,4),0.01*fs,0.005*fs,0.01*fs,fs);
figure;
surf(T2.*fs,F2,10*log10(P2),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of noisy signal from microphone 3');
%Spectrogram of beamformer output
[~,F3,T3,P3]=spectrogram(y_final,0.01*fs,0.005*fs,0.01*fs,fs);
figure;
surf(T3.*fs,F3,10*log10(P3),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of delay-and-sum beamformer output');

%% 2.1.A.3
%SNR of noisy signal of microphone 3
error_1=x(:,4)-source;
snr_1=snr(source,error_1);
fprintf('SNR of noisy signal = %4.2f dB\n',snr_1);
%SNR of beamformer output
error_2=y_final-source;
snr_2=snr(source,error_2);
fprintf('SNR of beamformer output = %4.2f dB\n',snr_2);

%% 2.1.B.1
source_cut=source(0.47*fs:0.5*fs);
x_3=x(0.47*fs:0.5*fs,4);
%Power spectrum of the frame from noisy signal
[P_x3,~]=pwelch(x_3,0.01*fs,0.005*fs,length(source_cut),fs,'twosided');
%Power spectrum of the frame from noise
[P_v3,f]=pwelch(x_3-source_cut,0.01*fs,0.005*fs,length(source_cut),fs,'twosided');
%Frequency response of Wiener filter
H_w=1-P_v3./P_x3;
%Plot of magnitude of the frequency response
figure;
plot(f/1000,10*log10(abs(H_w)));
xlim([0 8]);
xlabel( 'Frequency (kHz)');
ylabel ('Magnitude (dB)');
title ('Frequency response of IIR Wiener Filter');

%% 2.1.B.2
%Speech distortion index
nsd=(abs(1-H_w)).^2;
figure;
plot(f/1000,10*log10(nsd));
xlim([0 8]);
xlabel('Frequency (kHz)');
ylabel ('Magnitude (dB)');
title ('Speech Distortion Index');

%% 2.1.B.3
%Power spectrum of frame of clean signal
P_s=pwelch(source_cut,0.01*fs,0.005*fs,length(source_cut),fs,'twosided');
%Wiener filtering
S=fft(x_3);
Z=S.*H_w;
z=ifft(Z);
%Power spectrum of winer filter output
[P_z,f]=pwelch(z,0.01*fs,0.005*fs,length(source_cut),fs,'twosided');
%Plot of power spectra
figure;
plot(f/1000,10*log10(abs(P_s)));
hold on;
plot(f/1000,10*log10(abs(P_x3)));
hold on;
plot(f/1000,10*log10(abs(P_z)));
hold on;
plot(f/1000,10*log10(abs(P_v3)));
xlim ([0 8]);
legend ('Clear signal','Noisy signal','Output of Wiener filter','Noise');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Power Spectra'); 

%% 2.1.B.4
error_output=z-source_cut;
error_input=x_3-source_cut;
%SNR of wiener filter output
snr_3=snr(source_cut,error_input);
fprintf('SNR of wiener filter input = %4.2f dB\n',snr_3);
snr_4=snr(source_cut,error_output);
fprintf('SNR of wiener filter output = %4.2f dB\n',snr_4);
error_beamformer=y_final(0.47*fs:0.5*fs)-source_cut;
%SNR of beamformer output
snr_5=snr(source_cut,error_beamformer);
fprintf('SNR of beamformer output = %4.2f dB\n',snr_5);
%Power spectrum of frame of beamformer output
[P_b3,f]=pwelch(y_final(0.47*fs:0.5*fs),0.01*fs,0.005*fs,length(source_cut),fs,'twosided');
%Plot of power spectra
figure;
plot(f/1000,10*log10(abs(P_s)));
hold on;
plot(f/1000,10*log10(abs(P_x3)));
hold on;
plot(f/1000,10*log10(abs(P_z)));
hold on;
plot(f/1000,10*log10(abs(P_b3)));
xlim ([0 8]);
legend ('Clear signal','Noisy signal','Output of Wiener filter','Output of beamformer');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Power Spectra'); 

%% 2.2.A.1
N=7;
c=340;
d=0.04;
theta_s=pi/4;
fs=48000;
%Audioread of real signals
x_b(:,1)=audioread('sensor_0_2.wav');
x_b(:,2)=audioread('sensor_1_2.wav');
x_b(:,3)=audioread('sensor_2_2.wav');
x_b(:,4)=audioread('sensor_3_2.wav');
x_b(:,5)=audioread('sensor_4_2.wav');
x_b(:,6)=audioread('sensor_5_2.wav');
x_b(:,7)=audioread('sensor_6_2.wav');
source_b=audioread('source_2.wav');
X_b=fft(x_b);
X_b=fftshift(X_b,1);
d_ksb=zeros(length(x_b),7);
w =2*pi*fs*((-length(x_b)/2:length(x_b)/2-1)/length(x_b));
%Calculation of d(ks)
for i=0:N-1
    d_ksb(:,i+1)=exp(-1j*((N-1)/2)*(w/c)*d*cos(theta_s)).*exp(1j*i*(w/c)*d*cos(theta_s));
end
d_ksb=transpose(d_ksb);
d_ksb=ctranspose(d_ksb);
%Frequency response of beamformer
H_b=(1/N).*d_ksb;
Y_b=H_b.*X_b;
Y_b=ifftshift(Y_b,1);
y_b=ifft(Y_b);
%Beamformer output
y_finalb=sum(y_b,2);
y_finalb=real(y_finalb);
audiowrite('real_ds.wav',y_finalb,48000);
%% 2.2.A.2
t=(0:length(source_b)-1)*(1/fs);
%Plot of clean signal
figure;
plot(t,source_b);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Voice signal without noise");
%Plot of noisy signal
figure;
plot(t,x_b(:,4));
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Voice signal of microphone 3 with noise");
%Plot of beamformer output
figure;
plot(t,y_finalb);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Output of delay-and-sum beamformer");
%Spectrogram of clean signal
[~,F1b,T1b,P1b]=spectrogram(source_b,0.01*fs,0.005*fs,0.01*fs,fs);
figure;
surf(T1b.*fs,F1b,10*log10(P1b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of voice signal without noise');
%Spectrogram of noisy signal
[~,F2b,T2b,P2b]=spectrogram(x_b(:,4),0.01*fs,0.005*fs,0.01*fs,fs);
figure;
surf(T2b.*fs,F2b,10*log10(P2b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of noisy signal from microphone 3');
%Spectrogram of beamformer output
[~,F3b,T3b,P3b]=spectrogram(y_finalb,0.01*fs,0.005*fs,0.01*fs,fs);
figure;
surf(T3b.*fs,F3b,10*log10(P3b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of delay-and-sum beamformer output');

%% 2.2.A.3
%Frames of noisy signal
x_buf=buffer(x_b(:,4),0.03*fs);
[~,M]=size(x_buf);
%Noise
u=x_buf(:,15);
snr_w=zeros(1,M);
%SNR of each frame 
for i=1:M
    snr_w(i)=10*real(log10((mean(x_buf(:,i).^2)-mean(u.^2))/mean(u.^2)));
end
cnt=0;
for i=1:length(snr_w)
    if snr_w(i)>35
        snr_w(i)=35;
    elseif snr_w(i)<0
        snr_w(i)=0;
        cnt=cnt+1;
    end
end
%SSNR
ssnr=(1/(M-cnt)).*sum(snr_w);
fprintf('SSNR of noisy signal = %4.2f dB\n',ssnr);
%Frames of beamformer output
y_bufa=buffer(y_finalb,0.03*fs);
[~,M2]=size(y_bufa);
%Noise
u2=y_bufa(:,15);
snr_w2=zeros(1,M);
%SNR of each frame
for i=1:M
     snr_w2(i)=10*real(log10((mean(y_bufa(:,i).^2)-mean(u2.^2))/mean(u2.^2)));
end
cnt=0;
for i=1:length(snr_w2)
    if snr_w2(i)>35
        snr_w2(i)=35;
    elseif snr_w2(i)<0
        snr_w2(i)=0;
        cnt=cnt+1;
    end
end
%SSNR
ssnr2=(1/(M2-cnt)).*sum(snr_w2);
fprintf('SSNR of beamformer output = %4.2f dB\n',ssnr2);

%% 2.2.B.1
windowlength=0.03*fs;
window_overlap=0.012*fs;
windowstep=windowlength-window_overlap;
%Windowing the beamformer output
y_buff=buffer(y_finalb,windowlength,window_overlap);
y_buff=hamming(size(y_buff,1)).*y_buff;
%Power spectrum of noise
[P_ub,~]=pwelch(y_buff(:,15),0.01*fs,0.005*fs,0.03*fs,fs,'twosided');
%Power spectrum of filter input
[P_xb,~]=pwelch(y_buff,0.01*fs,0.005*fs,0.03*fs,fs,'twosided');
%Frequency response of wiener filter
H_wb=1-(P_ub./P_xb);
Y_buff=fft(y_buff);
Z_b=Y_buff.*H_wb;
%Filter output
z_b=ifft(Z_b);
[~,K]=size(z_b);
%Overlap-add method
z_finalb=z_b(window_overlap+1:end,1);
for i=1:K-1
    z_finalb=[z_finalb(1:end-window_overlap); z_b(windowstep+1:end,i)+z_b(1:window_overlap,i+1); z_b(window_overlap+1:end,i+1)];
end
z_finalb=z_finalb(1:length(source_b));
audiowrite('real_mmse.wav',z_finalb,48000);

%% 2.2.Â.2
t=(0:length(source_b)-1)*(1/fs);
%Plot of clean signal
figure;
plot(t,source_b);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Voice signal without noise");
%Plot of noisy signal
figure;
plot(t,x_b(:,4));
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Voice signal of microphone 3 with noise");
%Plot of wiener filter input
figure;
plot(t,y_finalb);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Input of Wiener filter");
ylim([-1 1]);
%Plot of wiener filter output
figure;
plot(t,z_finalb);
xlabel ("Time (sec)");
ylabel ("Amplitude");
title ("Output of Wiener filter");
ylim([-1 1]);
%Spectrogram of clean signal
[~,F1b,T1b,P1b]=spectrogram(source_b,0.03*fs,0.012*fs,0.03*fs,fs);
figure;
surf(T1b.*fs,F1b,10*log10(P1b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of voice signal without noise');
%Spectrogram of noisy signal
[~,F2b,T2b,P2b]=spectrogram(x_b(:,4),0.03*fs,0.012*fs,0.03*fs,fs);
figure;
surf(T2b.*fs,F2b,10*log10(P2b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of noisy signal from microphone 3');
%Spectrogram of wiener filter input
[~,F3b,T3b,P3b]=spectrogram(y_finalb,0.03*fs,0.012*fs,0.03*fs,fs);
figure;
surf(T3b.*fs,F3b,10*log10(P3b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of Wiener filter input');
%Spectrogram of wiener filter output
[~,F4b,T4b,P4b]=spectrogram(z_finalb,0.03*fs,0.012*fs,0.03*fs,fs);
figure;
surf(T4b.*fs,F4b,10*log10(P4b),'edgecolor','none');
axis tight; 
colormap(jet);
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram of Wiener filter output');

%% 2.2.B.3
%Frames of filter output
z_buf=buffer(z_finalb,0.03*fs);
[~,M3]=size(z_buf);
%Noise
u3=z_buf(:,15);
snr_w3=zeros(1,M3);
%SNR of each frame
for i=1:M3
    snr_w3(i)=10*real(log10((mean(z_buf(:,i).^2)-mean(u3.^2))/mean(u3.^2)));
end
cnt=0;
for i=1:length(snr_w3)
    if snr_w3(i)>35
        snr_w3(i)=35;
    elseif snr_w3(i)<0
        snr_w3(i)=0;
        cnt=cnt+1;
    end
end
%SSNR
ssnr3=(1/(M3-cnt)).*sum(snr_w3);
fprintf('SNR of Wiener filter output = %4.2f dB\n',ssnr3);

%% 2.2.B.4
ssnr_in=zeros(1,N);
%SSNR of noisy signal from each microphone
for i=1:N
    x_buf=buffer(x_b(:,i),0.03*fs);
    [~,M]=size(x_buf);
    u=x_buf(:,15);
    snr_in=zeros(1,M);
    for j=1:M
        snr_in(j)=10*real(log10((mean(x_buf(:,j).^2)-mean(u.^2))/mean(u.^2)));
    end
    cnt=0;
    for j=1:length(snr_in)
        if snr_in(j)>35
            snr_in(j)=35;
        elseif snr_in(j)<0
            snr_in(j)=0;
            cnt=cnt+1;
        end
    end
    ssnr_in(1,i)=(1/(M-cnt)).*sum(snr_in);  
end
%Mean of SSNRs
mean_ssnrin=mean(ssnr_in);
fprintf('Mean of SSNR of noisy signals = %4.2f dB\n',mean_ssnrin);
end