function [] = practice5() % Power Spectral Density
% - Generate a sumatorry of sine waves at 60,70 and 80 Hz with a sampling
% frequency of 1000 Hz for 0.1 seconds
% - Calculate the power spectral density using:
%   - periodogram (Non-Parametric)
%   - Yule Walker
%   - Burg
%   - Covariance 
%   - Modified covariance
% - Sketch the signal in time and frequency using FFT
% - Sketch the 5 graphs of PSD comparing the 5 methods
%       NOTE: Choose the order of the parametric methods.
% - Select a biological signal to obtain the PSD and observe the advantages

    % Original signal:
    Fs = 1000;
    t = 0:1/Fs:0.1;
    fn = [60,70,80];
    x = 0*t;
    for i = 1:length(fn), x = x + sin(2*pi*fn(i)*t); end
    % Fourier Transform:
    [f,X] = fourierVector(x,Fs);
    % Grafics:
    figure, subplot(2,1,1)
        plot(t,x)
        title('Sum of sine waves - Time domain')
        xlabel('Time (s)')
        ylabel('Amplitude')
    subplot(2,1,2)
        stem(f,X)
        title('Sum of sine waves - Frequency domain')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
    % Periodogram:
	[Pxx,Fp] = periodogram(x,[],Fs);
    Hp = dspdata.psd(Pxx,'Fs',Fs);
    figure, subplot(2,1,1)
        plot(Fp,Pxx)
        title('Periodogram')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
    subplot(2,1,2)
        plot(Hp)
    % Yule Walker:
    [YW,Fyw] = pyulear(x,26,[],Fs);
    Hyw = dspdata.psd(YW,'Fs',Fs);
    figure, subplot(2,1,1)
        plot(Fyw,YW)
        title('Yule Walker')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
    subplot(2,1,2)
        plot(Hyw)
    % Burg:
    [B,Fb] = pburg(x,26,[],Fs);
    Hb = dspdata.psd(B,'Fs',Fs);
    figure, subplot(2,1,1)
        plot(Fb,B)
        title('Burg')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
    subplot(2,1,2)
        plot(Hb)
    % Covariance:
    [C,Fc] = pcov(x,26,[],Fs);
    Hc = dspdata.psd(C,'Fs',Fs);
    figure, subplot(2,1,1)
        plot(Fc,C)
        title('Covariance')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
    subplot(2,1,2)
        plot(Hc)
    % Modified covariance:
    [MC,Fmc] = pmcov(x,26,[],Fs);
    Hmc = dspdata.psd(MC,'Fs',Fs);
    figure, subplot(2,1,1)
        plot(Fmc,MC)
        title('Modified covariance')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
    subplot(2,1,2)
        plot(Hmc)
        
	% Biological signal:
    [time,signal,Fs,~] = readPhysionet('ecg1');
    [fbio,Xbio] = fourierVector(signal,Fs);
    Mbio = pmcov(signal,26,[],Fs);
    Hbio = dspdata.psd(Mbio,'Fs',Fs);
    figure
        subplot(3,1,1)
            plot(time,signal)    
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('ECG signal')
            grid on
            xlim([0,3])
        subplot(3,1,2)
            plot(fbio,Xbio)
            title('Frequency domain')
            xlabel('Frequency (Hz)')
            ylabel('Amplitude')
        subplot(3,1,3)
            plot(Hbio)
            title('Power spectral density - Modified Covariance')
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
    X = abs(X(1:Ie))/length(X);
    f = (0:Ie-1)*((Fs/2)/(Ie-1));
end