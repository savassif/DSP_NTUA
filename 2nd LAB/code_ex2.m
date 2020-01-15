function code_ex2

%1.0

[music1,fs]=audioread('music-dsp19.wav'); 
%Calculation of the mean value of each row
music=sum(music1,2);
music=music./2;
%Normalization of music signal
normalized_music=music./max(abs(music));
%Time vector
t=(0:length(music)-1)*(1/fs);
%Plotting the normalized music signal
figure;
plot(t,normalized_music);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Normalized music signal');
R=2.^16;
%Partition of signal into sequences of 512 samples
sm_buf=buffer(normalized_music,512);

%1.1

function P=power_spectrum(signal)
   PN= 90.302;
   P= PN +10*log10((abs(fft(signal))).^2);
   P=P(1:length(P)/2);
end
%Adaptive quantizer
function q=quant(Bk, signal)
    range=max(signal)-min(signal);
    D=range./(2.^Bk);
    levels=zeros(1,2.^Bk);
    levels(1)=min(signal);
    q=zeros(1,length(signal));
    for i1=2:(2^Bk)
        levels(i1)=levels(i1-1)+D;
    end
    for i1=1:length(signal)
        diff=zeros(1,2.^Bk);
        for j1=1:length(levels)
            diff(j1)=abs(signal(i1)-levels(j1));
        end
        q(i1)=levels(diff==min(diff));
    end
end    
%Non-Adaptive quantizer
function q_st=quant_st(signal)
    range=1-(-1);
    Bits=8;
    D=range/(2^Bits);
    levels=zeros(1,2^Bits);
    levels(1)=-1;
    q_st=zeros(1,length(signal));
    for i1=2:(2^Bits)
        levels(i1)=levels(i1-1)+D;
    end
    for i1=1:length(signal)
        diff=zeros(1,2^Bits);
        for j1=1:length(levels)
            diff(j1)=abs(signal(i1)-levels(j1));
        end
        q_st(i1)=levels(diff==min(diff));
    end
end

%Frequency vector
f1=fs*(1:256)/512;
%Absolute threshold of hearing
Tq= @(f)3.64*((f./1000).^(-0.8)) -6.5*exp(-0.6*(f./1000-3.3).^2) + (10.^(-3))*((f./1000).^4);
Tq_0 = Tq(f1);
%Bark scale
b= @(f) 13*atan(0.00076.*f) + 3.5*atan((f./7500).^2);
b_0=b(f1);
signal_final_st=0;

function SF=SF(delta, p)
    if (delta>=-3)&&(delta<-1)
        sf=17*delta-0.4.*p+11;
    else
        if (delta>=-1)&&(delta<0)
        sf=(0.4*p+6)*delta;
        else
            if (delta>=0)&&(delta<1)
        sf=-17*delta; 
        else
            sf=(0.15*p-17)*delta-0.15*p;
            end
        end
    end
    SF=sf;
end

M=32;
N=2*M;
n=(0:N-1);
h=zeros(N,M);
g=zeros(size(h));
center_freq=zeros(1,M);
signal_final=zeros(650240,1);
%Filters
for k=1:M
    h(:,k)=(sin((n+1/2).*(pi/(2*M)))).*(sqrt(2/M)).*cos(((2*n+M+1).*(2*(k-1)+1)*pi)./(4*M));
    g(:,k)=(sin((N-1-n+1/2).*(pi/(2*M)))).*(sqrt(2/M)).*cos(((2*(N-1-n)+M+1).*(2*(k-1)+1)*pi)./(4*M));
    center_freq(1,k)= (2*k-1)*fs * pi/(2*M);
end
%Plots of filter impulse responses for k=10
figure;
plot (n,h(:,10));
xlabel('n');
ylabel('Amplitude');
title('Filter impulse response');
figure;
plot(n,g(:,10));
xlabel('n');
ylabel('Amplitude');
title('Filter impulse response');
for i=1:size(sm_buf,2)
    %Windowing the signal with hanning window
    windowed_signal=hanning(512).*sm_buf(:,i);
    %Power spectrum
    P1=power_spectrum(windowed_signal);
    Ptm=zeros(256,1);
    St=zeros(256,1);
    Tg=zeros(256,1);
    %Tone maskers 
    [maximum, points]=findpeaks(P1(3:250));
    for j=1:length(points)
        maxk=maximum(j);
        k=points(j)+2; %match P1 indices
        if ((P1(k-2)+7)<maxk)&&((P1(k+2)+7)<maxk)
            if k<63
                St(k)=1;
                Ptm(k)=10*log10(10.^(0.1*P1(k-1))+10.^(0.1*P1(k))+10.^(0.1*P1(k+1)));
            else
                if ((P1(k-3)+7)<maxk)&&((P1(k+3)+7)<maxk)
             if (k>=63)&&(k<127)
                    St(k)=1;
                    Ptm(k)=10*log10(10.^(0.1*P1(k-1))+10.^(0.1*P1(k))+10.^(0.1*P1(k+1)));
                else
                    if ((P1(k-4)+7)<maxk)&&((P1(k+4)+7)<maxk)&&((P1(k-5)+7)<maxk)&&((P1(k+5)+7)<maxk)&&((P1(k-6)+7)<maxk)&&((P1(k+6)+7)<maxk)
                        St(k)=1;  
                        Ptm(k)=10*log10(10.^(0.1*P1(k-1))+10.^(0.1*P1(k))+10.^(0.1*P1(k+1)));
                    end
             end
                end
            end
        end
    end
    %Noise maskers
    Pnm=findNoiseMaskers(P1, Ptm, b_0);
    %Plots for the 100th sequence
    if i==100
        %Plot of tone maskers along with power spectrum
        Ptmplot=Ptm;
        Ptmplot(Ptmplot==0)=NaN;
        figure;
        plot(b_0, P1);
        hold on;
        stem(b_0, Ptmplot);
        xlabel('Frequency (Bark)');
        ylabel ('Magnitude');
        title('Tone maskers');
        %Plot of noise maskers along with power spectrum
        Pnmplot=Pnm;
        Pnmplot(Pnmplot==0)=NaN;
        figure;
        plot(b_0, P1);
        hold on;
        stem(b_0, Pnmplot);
        xlabel('Frequency (Bark)');
        ylabel ('Magnitude');
        title('Noise maskers');
    end
    %Checked Tone and Noise maskers
    [Ptm_new, Pnm_new]=checkMaskers(Ptm',Pnm,Tq_0,b_0);
    %Plots for the 100th sequence 
    if i==100
        %Plot of checked tone maskers along with power spectrum
        Ptm_newplot=Ptm_new;
        Ptm_newplot(Ptm_newplot==0)=NaN;
        figure;
        plot(b_0, P1);
        hold on;
        stem(b_0, Ptm_newplot);
        xlabel('Frequency (Bark)');
        ylabel ('Magnitude');
        title('Checked tone maskers');
        %Plot of checked noise maskers along with power spectrum
        Pnm_newplot=Pnm_new;
        Pnm_newplot(Pnm_newplot==0)=NaN;
        figure;
        plot(b_0, P1);
        hold on;
        stem(b_0, Pnm_newplot);
        xlabel('Frequency (Bark)');
        ylabel ('Magnitude');
        title('Checked noise maskers');
    end
    %Indices of tone and noise maskers
    indices_nm=find(Pnm_new);
    indices_tm=find(Ptm_new);
    Sf_nm=zeros(256,256);
    Sf_tm=zeros(256,256);
    T_nm=zeros(size(Sf_nm));
    T_tm=zeros(size(Sf_tm));
    for cntj=1:length(indices_nm)
        %Indices of 12 Bark region 
       indices_b_nm=find((b_0>=b_0(indices_nm(cntj))-3)& (b_0<=b_0(indices_nm(cntj))+8));
       for cnti=1:length(indices_b_nm)
          deltab_nm=b_0(indices_b_nm(cnti))-b_0(indices_nm(cntj));
          Sf_nm(indices_b_nm(cnti),indices_nm(cntj))=SF(deltab_nm,Pnm_new(indices_nm(cntj)));
          %Noise masking thresholds
          T_nm(indices_b_nm(cnti),indices_nm(cntj))=Pnm_new(indices_nm(cntj))-0.175*b_0(indices_nm(cntj))+Sf_nm(indices_b_nm(cnti),indices_nm(cntj))-2.025;
       end 
    end
    for cntj=1:length(indices_tm)
        %Indices of 12 Bark region  
       indices_b_tm=find((b_0>=b_0(indices_tm(cntj))-3)& (b_0<=b_0(indices_tm(cntj))+8));
       for cnti=1:length(indices_b_tm)
          deltab_tm=b_0(indices_b_tm(cnti))-b_0(indices_tm(cntj));
          Sf_tm(indices_b_tm(cnti),indices_tm(cntj))=SF(deltab_tm,Ptm_new(indices_tm(cntj)));
          %Tone masking thresholds
          T_tm(indices_b_tm(cnti),indices_tm(cntj))=Ptm_new(indices_tm(cntj))-0.275*b_0(indices_tm(cntj))+Sf_tm(indices_b_tm(cnti),indices_tm(cntj))-6.025;
       end  
    end
    %Plots of tone and noise masking thresholds
    if i==100
        figure;
        surf(b_0,b_0 , T_tm, 'edgecolor', 'none','FaceAlpha',0.5);
        xlabel('Frequency (Bark)(j)');
        ylabel ('Frequency (Bark)(i)');
        title('Tone masking thresholds');
        figure;
        surf(b_0,b_0, T_nm,'edgecolor','none','FaceAlpha',0.5);
        xlabel('Frequency (Bark)(j)');
        ylabel ('Frequency (Bark)(i)');
        title('Noise masking thresholds');
    end
    %Sums calculated for global masking threshold
    for cnt1=1:256
        sum_1=0;
        sum_2=0;
        for cnt3=1:256
            if T_tm(cnt1,cnt3)~=0
                sum_1=sum_1+10.^(0.1.*(T_tm(cnt1,cnt3)));
            end
            if T_nm(cnt1,cnt3)~=0
                sum_2=sum_2+10.^(0.1.*T_nm(cnt1,cnt3));
            end
        end
        %Global Masking Threshold
        Tg(cnt1,1)=10.*log10(10.^(0.1.*Tq_0(cnt1))+sum_1+sum_2);
    end
    %Plot of absolute threshold of hearing along with global masking
    %threshold for the 100th sequence
    if i==100
        figure;
        plot (b_0, Tq_0);
        hold on;
        plot(b_0, Tg);
        xlabel('Frequency (Bark)');
        ylabel('Magnitude');
        title ('Thresholds');
        legend ('Absolute Threshold of Hearing','Global Masking Threshold');
    end
    %2nd part
    u=zeros(512+N-1,M);
    y=zeros((512+N)/M,M);
    B=zeros(1,M);
    xq=zeros(18,M);
    xq_st=zeros(18,M);
    w=zeros(size(y,1)*M,M);
    w_st=zeros(size(w));
    x_final=zeros(size(w,1)+N-1,M);
    x_st_final=zeros(size(x_final,1));
    for k=1:M
        %Convolution of signal with impulse response h
        u(:,k)=conv(sm_buf(:,i),h(:,k));
        %Downsampling by factor 32
        y(:,k)=downsample(u(:,k),M);
        %Indices for the right choice of Tg values
        indices_1=(8*(k-1)+1):8*k;
        %Number of quantization bits
        B(1,k)=ceil(log2(R/min(Tg(indices_1)))-1);
        %Quantized signal with adaptive quantizer
        xq(:,k)=quant(B(1,k),y(:,k));
        %Quantized signal with non-adaptive quantizer
        xq_st(:,k)=quant_st(y(:,k));
        %Upsampling by factor 32
        w(:,k)=upsample(xq(:,k),M);
        w_st(:,k)=upsample(xq_st(:,k),M);
        %Convolution of w with impulse response g
        x_final(:,k)=conv(w(:,k),g(:,k));
        x_st_final(:,k)=conv(w_st(:,k),g(:,k));
    end
    %Reconstructed sequence
    s=sum(x_final,2);
    s_st=sum(x_st_final,2);
    %Overlap-Add Method
    if i==1
        signal_final=s;
        signal_final_st=s_st;
    else
        signal_final=[signal_final(1:end-2*N+1); signal_final(end-2*N+2:end)+s(1:end-512); s(end-512+1:end)];
        signal_final_st=[signal_final_st(1:end-2*N+1); signal_final(end-2*N+2:end)+s(1:end-512); s(end-512+1:end)];
    end
    %Compressio ratio of adaptive quantizer for the 100th sequence
    if i==100
        comp_ratio=(mean(B)./16)*100;
    end
end
%Shifting the signal by N=64 samples
signal_final=signal_final(N:length(music)+N-1);
signal_final_st=signal_final_st(N:length(music)+N-1);
%Compression ratio of non-adaptive quantizer
comp_ratio_st=(8/16)*100;
%Error (adaptive quantizer)
error=normalized_music-signal_final;
%Error (non-adaptive quantizer)
error_st=normalized_music-signal_final_st;
%Mean Squared Error (adaptive quantizer)
mse=immse(normalized_music,signal_final);
%Mean Squared Error (non-adaptive quantizer)
mse_st=immse(normalized_music,signal_final_st);
%Plot of reconstructed music signal after adaptive quantising
figure;
plot(t,signal_final);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Final Signal after Adaptive Qunatising');
%Plot of reconstructed music signal after non-adaptive quantising
figure;
plot(t,signal_final_st);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Final Signal after Non-Adaptive Qunatising');
%Plot of error of adaptive quantising
figure;
plot(t,error);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Error of Adaptive Quantising');
%Plot of error of non-adaptive quantising
figure;
plot(t,error_st);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Error of Non-Adaptive Quantising');
audiowrite('music_final.wav',signal_final,fs);
audiowrite('music_final_st.wav',signal_final_st,fs);
end
