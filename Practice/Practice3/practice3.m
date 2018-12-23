function practice3() % Digital filters
%% Part I: Design of a band-pass filter from an analog filter

% Obtain the bode diagrams of the analog and digital filters.   
    Fs = 1500;
    n = 4096; % Number of points to be ploted
    figure
        freqz([0.3249,0.3249],[1.3249,-0.6751],n,Fs)
        title('Low pass 150Hz filter')
    figure
        freqz([1,-1],[1.0021,-0.9979],n,Fs)
        title('High pass 1Hz filter')
    B = [0.3249,0,-0.3249];
    A = [1.3277,-1.9986,0.6737];
    figure
        freqz(B,A,n,Fs)
        title('Band pass 1-150Hz filter')

% Select an ECG from PhysioNet and apply your filter to it:
    [time,signal,~,~] = readPhysionet('ecg1');
    fsignal = filter(B,A,signal);
    figure
        subplot(2,1,1)
            plot(time,signal)
            title('Original ECG signal')
            xlabel('Time (s)')
            ylabel('Amplitude')
            xlim([0,2])
            grid on
        subplot(2,1,2)
            plot(time,fsignal)
            title('Filtered ECG signal')
            xlabel('Time (s)')
            ylabel('Amplitude')
            xlim([0,2])
            grid on

% Test your filter with a summatory of sine-waves to observe its behaviour 
% in the time and frequency domain.
    F = [1,10,50,100,150,200];
    [t,signal2] = sumSines(F,Fs);
    fsignal2 = filter(B,A,signal2);
    [f2,S2] = fourierVector(signal2,Fs);
    [ff2,FS2] = fourierVector(fsignal2,Fs);
    figure
        subplot(2,2,1)
            plot(t,signal2)
            title('Original sum of sines signal')
            xlabel('Time (s)')
            ylabel('Amplitude')
            xlim([0,0.2])
            grid on
        subplot(2,2,2)
            plot(t,fsignal2)
            title('Filtered sum of sines signal')
            xlabel('Time (s)')
            ylabel('Amplitude')
            xlim([0,0.02])
            grid on
        subplot(2,2,3)
            plot(f2,S2)
            title('Frequency domain original signal')
            xlabel('Frequency (Hz)')
            ylabel('Amplitude')
            xlim([0,250])
            grid on
        subplot(2,2,4)
            plot(ff2,FS2)
            title('Frequency domain filtered signal')
            xlabel('Frequency (Hz)')
            ylabel('Amplitude')
            xlim([0,250])
            grid on
    
%% Part II: IIR and FIR filters

% Design a low pass filter in the 4 diferent configurations (Butterwooth,
% Cheby I, Cheby II, Elliptic) fc = 60 Hz, order 3.
    fc = 60;
    fs = 600;
    [b,a] = butter(3,fc/(fs/2));
    figure
        freqz(b,a,n,fs)
        title('Butter 60 Hz low pass filter')
    [b,a] = cheby1(3,3,fc/(fs/2));
    figure
        freqz(b,a,n,fs)
        title('Cheby1 60 Hz low pass filter 3 dB Bandpass Ripple')
    [b,a] = cheby2(3,30,fc/(fs/2));
    figure
        freqz(b,a,n,fs)
        title('Cheby2 60 Hz low pass filter 30 dB StopBand Ripple')
    [b,a] = ellip(3,3,30,fc/(fs/2));
    figure
        freqz(b,a,n,fs)
        title('Elliptic 60 Hz low pass 3db Bandpass and 30dB Stopband Ripple')
    
% Choose 5 different FIR filters and compare their behaviour. Find the
% order to show a similar behaviour of an IIR filter.
    [b,a] = fir1(6,fc/(fs/2),hamming(7));
    figure
        freqz(b,a,n,fs)
        title('FIR filter Hamming window')
	[b,a] = fir1(6,fc/(fs/2),hann(7));
    figure
        freqz(b,a,n,fs)
        title('FIR filter Hann window')
    [b,a] = fir1(6,fc/(fs/2),gausswin(7));
    figure
        freqz(b,a,n,fs)
        title('FIR filter Gauss window')
    [b,a] = fir1(6,fc/(fs/2),kaiser(7));
    figure
        freqz(b,a,n,fs)
        title('FIR filter Kaiser window')
	[b,a] = fir1(8,fc/(fs/2),blackmanharris(9));
    figure
        freqz(b,a,n,fs)
        title('FIR filter Blackman-Harris window')
    
% Select a biosignal from PhysioNet and design a digital filter to
% observe a specifdic characteristic of your signal. Explain why did you
% choose that filter
	[b,a] = fir1(15,30/(1000/2),blackmanharris(16));
        figure
        freqz(b,a,4095,1000)
    title('FIR filter Blackman-Harris window for abdominal (fetus) ECG')    
    [time,signals,~,~] = readPhysionet('f_a_ecg');
    signal3 = signals(:,2);
    fsignal3 = filter(b,a,signal3);
    figure
    subplot(2,1,1)
        plot(time,signal3)
        title('Original abdominal (fetus) ECG signal')
        xlabel('Time (s)')
        ylabel('Amplitude')
        xlim([0,1])
        ylim([-70,30])
        grid on
    subplot(2,1,2)
        plot(time,fsignal3)
        title('Filtered abdominal (fetus) ECG signal')
        xlabel('Time (s)')
        ylabel('Amplitude')
        xlim([0,1])
        ylim([-70,30])
        grid on

end

function [time,val,Fs,labels] = readPhysionet(Name)
% This function reads a pair of files (Name.mat and Name.info) from a 
% PhysioBank record, baseline-corrects and scales the time series contained
% in the .mat file.  The baseline-corrected and scaled time series are the 
% rows of matrix 'val', and each column contains simultaneous samples of 
% each time series.

% Read mat File:
    matName = strcat(Name, '.mat');
    load(matName);
    n = size(val,1);
% Read info File:
    infoName = strcat(Name, '.info');
    fid = fopen(infoName, 'rt');
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    freqint = sscanf(fgetl(fid), 'Sampling frequency: %f Hz  Sampling interval: %f sec');
    Fs = freqint(1);
    interval = freqint(2);
    fgetl(fid);
% Read data of each signal
    signal = cell(1,n);
    gain = zeros(1,n);
    base = zeros(1,n);
    units = cell(1,n);
    for i = 1:n
      [~, signal(i), gain(i), base(i), units(i)] = strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
    fclose(fid);

% Baseline-corrects and scales the time series:
    val(val==-32768) = NaN;
    for i = 1:n
        val(i, :) = (val(i, :) - base(i)) / gain(i);
    end
    time = (1:size(val, 2)) * interval;
    val = val';
    
% Gives information of each signal:
    labels = cell(1,length(signal));
    for i = 1:length(signal)
        labels{i} = strcat(signal{i}, ' (', units{i}, ')'); 
    end
    
end

function [f,X] = fourierVector(x,Fs)
% To obtain absolute positive values of the fourier transform of x
% it also obtains the corresponding frequency vector
    X = fft(x);
    Ie = floor(length(X)/2);
    X = abs(X(:,1:Ie))/length(X);
    f = (0:Ie-1)*((Fs/2)/(Ie-1));
end

function [t,x] = sumSines(F,Fs)
    t = 0:1/Fs:1-1/Fs;
    x = 0*t;
    for i = 1:length(F)
        x = x + sin(2*pi*F(i)*t);
    end
end