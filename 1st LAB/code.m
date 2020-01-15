function code
%PART 1    
    function indices=mypeaks(x,tones_numb)
        for cnt=1:tones_numb
        [pks,locs]=findpeaks(x(cnt,:));
         [~, I] = maxk(pks,4);
         I=sort(I);
         ind(cnt,:)= locs(I(3:4));
        end
        indices = ind;
    end

%1.1
%Create Tones
n = linspace(1,1000,1000);
d0=sin(0.7217.*n) + sin(1.0247.*n);
d1=sin(0.5346.*n) + sin(0.9273.*n);
d2=sin(0.5346.*n) + sin(1.0247.*n);
d3=sin(0.5346.*n) + sin(1.1328.*n);
d4=sin(0.5906.*n) + sin(0.9273.*n);
d5=sin(0.5906.*n) + sin(1.0247.*n);
d6=sin(0.5906.*n) + sin(1.1328.*n);
d7=sin(0.6535.*n) + sin(0.9273.*n);
d8=sin(0.6535.*n) + sin(1.0247.*n);
d9=sin(0.6535.*n) + sin(1.1328.*n);

%1.2
% W Vector
w =2*pi*((-length(d4)/2:length(d4)/2-1)/length(d4));

%DFT of d4 and d6
D4=fft(d4);
D4=fftshift(D4)/1000;
plot(w,abs(D4));
xlabel('Frequency');
ylabel('Amplitude');
title('DFT of d4');

figure;

D6=fft(d6);
D6=fftshift(D6)/1000;
plot(w,abs(D6));
xlabel('Frequency');
ylabel('Amplitude');
title('DFT of d6');

%1.3
zero_vector=zeros(1,100);

%adding zeros to Signal
toneseq=[d0 zero_vector d6 zero_vector d2 zero_vector d3 zero_vector d2 zero_vector d1 zero_vector d0  zero_vector d6];
audiowrite('tone_sequence.wav',toneseq/max(abs(toneseq)),8192);

%Creation of Windows and Windowed Signal(s)
rect=rectwin(1000);
j=1;
windowedSignal1 = zeros(8,1000);

for i=1:8
    windowedSignal1(i,:)=rect'.*toneseq(1, j:(j-1+1000));
    j=j+1100;
end

ham=hamming(1000);
j=1;
windowedSignal2 = zeros(8,1000);

for i=1:8
    windowedSignal2(i,:)=ham'.*toneseq(1, j:(j-1+1000));
    j=j+1100;
end

WindowedSignal1 = zeros(8,1000);
WindowedSignal2 = zeros(8,1000);

%DFT of Windowed Signals
for i=1:8
WindowedSignal1(i,:)=fft(windowedSignal1(i,:));
WindowedSignal2(i,:)=fft(windowedSignal2(i,:));
WindowedSignal1(i,:)=fftshift(WindowedSignal1(i,:))/1000;
WindowedSignal2(i,:)=fftshift(WindowedSignal2(i,:))/1000;
end

%1.5 
%List of Indices k

all_tones=[d0; d1; d2; d3; d4; d5; d6; d7; d8; d9];

%DFT of all_tones
All_tones=zeros(size(all_tones));
for i=1:10
    All_tones(i,:)=fft(all_tones(i,:));
    All_tones(i,:)=fftshift(All_tones(i,:))/1000;
end

%list of indices k and corresponding frequencies
points = mypeaks(abs(All_tones),10);
frequencies=w(points);
k=[points(:,1)', points(:,2)'];

for i=1:10
    points(i,3)=i-1;
end

k=unique(k);
k(2,:)=w(k(1,:));

%1.6
function vec=ttdecode(x)
   l=1;
   %We find the tones of the signal and store them in the array "tones"
   while (i<=length(x)) 
        if x(i)~=0
            tones(l,:)=x(i:i+999);
            l=l+1;
            i=i+1000;
        else
            i=i+1;
        end
   end
   [numb_tones , ~]=size(tones);
   Tones=zeros(size(tones));
   %DFT of the tones
   for i=1:numb_tones
      Tones(i,:)=fft(tones(i,:));
      Tones(i,:)=fftshift(Tones(i,:))/1000;
   end 
   %Energy calculation of each tone
   Energy=(abs(Tones)).^2;
   locs=mypeaks(Energy,numb_tones);
   vector=zeros(1,numb_tones);
   for j=1:numb_tones
       for i=1:10
           if locs(j,:)==points(i,1:2)
               vector(1,j)=points(i,3);
               break;
           end
       end
   end
   vec=vector;
   disp(vector);
end

number=ttdecode(toneseq);

%1.7
S=load('my_touchtones.mat');
number1=ttdecode(S.easySig);
number2=ttdecode(S.hardSig);

%PART 2

%2.1

Fs=1000;
n=0:(1/Fs):2;
%White Gaussian noise
v=randn(1,length(n));
%1st Signal
x=2*cos(2*pi*70*n) + 3*sin(2*pi*140*n) + 0.15*v;
figure;
plot(n,x);
xlabel('Time(s)');
ylabel('Amplitude');
title('1st Signal');

figure;

%STFT
[S1,F1,T1]=spectrogram(x, 0.04*Fs, 20,40,Fs);
surf(T1.*1000,F1,abs(S1),'edgecolor','none');
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram using STFT');

figure;
%Scale and Frequency Calculation
[s1,f1]=wavescales('morl',Fs);
%DT-CWT
cwtstruct1 = cwtft({x,1/Fs},'Scales',s1,'Wavelet','morl');
cfs1 = cwtstruct1.cfs;
surf(n.*1000,f1',abs(cfs1),'edgecolor','none');
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram using Wavelets');

%2.2

%White Gaussian noise
v=randn(1,length(n));
%2nd Signal
y=1.7*cos(2*pi*90*n) + 0.15*v+1.7*double(eq(n,0.625))+1.7*double(eq(n,0.800));
figure;
plot (n,y);xlabel('Time(s)');
ylabel('Amplitude');
title('2nd Signal');

%STFT
[S2,F2,T2]=spectrogram(y, 0.04*Fs, 20,40,Fs);
figure;
contour(T2.*1000,F2,abs(S2));
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram using STFT');

%Scale and Frequency Calculation
[s2,f2]=wavescales('morl',Fs);
cwtstruct2 = cwtft({y,1/Fs},'Scales',s2,'Wavelet','morl');
cfs2 = cwtstruct2.cfs;
figure;
contour(n.*1000,f2',abs(cfs2));
view(0,90);
xlabel('Time(ms)');
ylabel('Frequency(Hz)');
title('Spectrogram using Wavelets');

%PART 3

[speechSignal, Fs1] = audioread('speech_utterance.wav');
[musicSignal, Fs2]=audioread('music_cut.wav');

    function ste=STE (signal, Fs, timewinlength)
        winLen = timewinlength*Fs;
        winOverlap = winLen-1;
        % A hamming window is chosen
        wHamm = hamming(winLen);
        
        % Framing and windowing the signal
        sigFramed = buffer(signal, ceil(winLen), ceil(winOverlap), 'nodelay');
        sigWindowed = diag(sparse(wHamm)) * sigFramed;

        % Short-Time Energy calculation
        energyST = sum(sigWindowed.^2,1);
        ste=energyST;
    end

    function zcr=ZCR(signal, Fs, timewinlength)
        winLen = timewinlength*Fs;
        winOverlap = winLen-1;
        % A hamming window is chosen
        wHamm = hamming(winLen);
        
        signs=abs(diff(signal>=0));
        % Framing and windowing the signal "signs"
        signsFramed = buffer(signs, ceil(winLen), ceil(winOverlap), 'nodelay');
        signsWindowed = diag(sparse(wHamm)) * signsFramed;
        % Zero Crossing Rate calculation
        ZCR= mean(signsWindowed,1);
        zcr=ZCR;
    end

    function common_plots(func1,func2,Fs,delay)
        t=(0:length(func1)-1)/Fs;
        figure;
        %For the speech signal, plots are not normalized
        if Fs==Fs1
            plot(t, func1);
            hold on;   
            plot(t(floor(delay)+1:end - round(delay)), func2, 'r');
        %For the music signal, plots are normalized 
        else 
            plot(t, func1/max(abs(func1)));
            hold on;   
            plot(t(floor(delay)+1:end - round(delay)), func2/max(abs(func2)), 'r');
        end
        xlabel('Time (sec)');
    end

    function procedure_ste(signal,Fs,timewindowlength,type)
        %Short-Time Energy calculation and common plot of signal and short-time
        %energy 
        winlen=timewindowlength*Fs; 
        % Short-Time energy is delayed due to lowpass filtering. This delay is
        % compensated for the graph.
        delay_ste = (winlen-1)/2;
        short_time_energy=STE(signal, Fs, timewindowlength);
        common_plots(signal,short_time_energy,Fs, delay_ste);
        legend({type,'Short-Time Energy'});
        title([type, ' Signal Along With Short-Time Energy with window length ', num2str(timewindowlength),' sec']);
    end

    function multiple_graphs_ste(signal,Fs,windowlength,type)
        % Short-Time Energy calculation and plots for various window length
        % values
        figure;
        t = (0:length(signal)-1)/Fs;
        for i=1:length(windowlength)
            short_time_energy=STE(signal, Fs, windowlength(i));
            winLen = ceil(windowlength(i)*Fs);
            delay = (winLen-1)/2;
            subplot(length(windowlength),1,i);
            plot(t(floor(delay)+1:end - round(delay)), short_time_energy);
            xlim([0 t(end)]);
            xlabel('Time (sec)');
            title(['STE of ', type,' Signal with window length ', num2str(windowlength(i)),'sec']);
        end
    end

    function procedure_zcr(signal, Fs, timewindowlength,type)
        winlen=timewindowlength*Fs; 
        % Zero Crossing Rate is delayed due to lowpass filtering. This delay is
        % compensated for the graph.
        delay_zcr = (winlen)/2;
        zero_crossing_rate=ZCR(signal, Fs, timewindowlength);
        common_plots(signal,zero_crossing_rate,Fs, delay_zcr);
        legend({type,'Zero Crossing Rate'});
        title([type, ' Signal Along With Zero Crossing Rate with window length ', num2str(timewindowlength),' sec']);
       
    end

    function multiple_graphs_zcr(signal,Fs,windowlength,type)
        % Zero Crossing Rate calculation and plots for various window length
        % values
        figure;
        t = (0:length(signal)-1)/Fs;
        for i=1:length(windowlength)
           zero_crossing_rate=ZCR(signal, Fs, windowlength(i));
            winLen = ceil(windowlength(i)*Fs);
            delay = winLen/2;
            subplot(length(windowlength),1,i);
            plot(t(floor(delay)+1:end - round(delay)), zero_crossing_rate);
            xlim([0 t(end)]);
            xlabel('Time (sec)');
            title(['ZCR of ', type, ' Signal with window length ', num2str(windowlength(i)),'sec'])
        end
    end

%Short-Time Energy calculation and plots for speech signal with window length 20 ms
procedure_ste(speechSignal,Fs1,0.02,'Speech');
%Short-Time Energy calculation and plots for speech signal with window length 30 ms
procedure_ste(speechSignal,Fs1,0.03,'Speech');

%Zero Crossing Rate calculation for speech signal with window length 20 ms
procedure_zcr(speechSignal,Fs1,0.02,'Speech');
%Zero Crossing Rate calculation for speech signal with window length 30 ms
procedure_zcr(speechSignal,Fs1,0.03,'Speech');

windowlengths=[0.01 0.026 0.06 0.1];
%Speech signal
%Calculation and plots of short-time energy for various window lengths
multiple_graphs_ste(speechSignal, Fs1, windowlengths,'Speech');

%Calculation and plots of zero crossing rate for various window lengths
multiple_graphs_zcr(speechSignal, Fs1, windowlengths, 'Speech');


%Short-Time Energy calculation and plots for music signal with window length 20 ms
procedure_ste(musicSignal,Fs2,0.02,'Music');
%Short-Time Energy calculation and plots for music signal with window length 30 ms
procedure_ste(musicSignal,Fs2,0.03,'Music');

%Zero Crossing Rate calculation for music signal with window length 20 ms
procedure_zcr(musicSignal,Fs2,0.02,'Music');
%Zero Crossing Rate calculation for music signal with window length 30 ms
procedure_zcr(musicSignal,Fs2,0.03,'Music');

%Music signal
%Calculation and plots of short-time energy for various window lengths
multiple_graphs_ste(musicSignal, Fs2, windowlengths,'Music');

%Calculation and plots of zero crossing rate for various window lengths
multiple_graphs_zcr(musicSignal, Fs2, windowlengths,'Music');
end 
