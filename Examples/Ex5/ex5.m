function [] = ex5()

% Read and plot signal:
    [time,signal,labels] = readPhysionet('f_a_ecg');
    figure, plot(time,signal); legend(labels); xlabel('Time(s)');
    
% Filters:
    fs = 100;
    [B,A] = butter(3,[2,20]/(fs/2)); % Coeficients
    fsignal = filter(B,A,signal); % Aply filter
    ffsignal = filtfilt(B,A,signal); % Aply doble filter
    
% Fourier Transform:
    [X1,freqs] = fourierVector(signal,fs);
    [X2,freqf] = fourierVector(fsignal,fs);
    [X3,freqff] = fourierVector(ffsignal,fs);

% Graphics:
    for i = 1:size(signal,2)
        figure 
            subplot(2,1,1)
                plot(time,signal(:,i),'r')
                hold on
                plot(time,fsignal(:,i),'b')
                plot(time,ffsignal(:,i),'g')
                title(['Time domain ',labels{i}])
                xlabel('Time (s)')
                legend('Signal','filt signal','filtfilt signal')
                xlim([0,1])
            subplot(2,1,2)
                plot(X1,freqs(:,i),'r')
                hold on
                plot(X2,freqf(:,i),'b')
                plot(X3,freqff(:,i),'g')
                title(['Frequency domain ',labels{i}])
                xlabel('Frequency (Hz)')
                legend('Signal','filt signal','filtfilt signal')
    end

end

function [time,val,labels] = readPhysionet(Name)

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
    X = abs(X(1:Ie,:))/length(X);
    f = (0:Ie-1)*((Fs/2)/(Ie-1));
end
