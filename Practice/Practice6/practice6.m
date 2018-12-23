function [] = practice6()
    % Signal:
    t = 0:0.001:2;
    x = chirp(t,0,1,150);

	% Biological signal:
    [time,signal,Fs,~] = readPhysionet('ecg1');
        
	% Wavelet parameters:
%     M = 1; % Number of voices
%     a = [];
%     a1 = zeros(1,M);
%     for j = 1:20 % Generate Scales:
%         for m = 1:M
%             a1(1,m) = 2^(j+(m/M));
%         end
%         a = [a,a1];
%     end
    a = 1:100;
    
    % Wavelet for signal:
    c = cwt(hilbert(x),a,'db3'); % Change mother wavelet
    % c = cwt(x,a,'db3'); % Change mother wavelet
    f = ones(1,length(a))*1000./a; % Convert scale-frequency
    
	% Wavelet for biosignal:
    bc = cwt(hilbert(signal),a,'haar'); % Change mother wavelet
    % c = cwt(x,a,'db3'); % Change mother wavelet
    bf = ones(1,length(a))*Fs./a; % Convert scale-frequency
    
	% Graphics:
    figure,
        subplot(2,2,1)
            plot(t,x)
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('Chirp signal')
        subplot(2,2,2)
            plot(time,signal)    
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('ECG signal')
            xlim([0,5])
    	subplot(2,2,3)
            imagesc(t,f,abs(c));
            colormap(jet);
            axis('ij');
            xlabel('Time (s)')
            ylabel('Frequency')
            title('Chirp signal')
        subplot(2,2,4)
            imagesc(time,bf,abs(bc));
            colormap(jet);
            axis('ij');
            xlabel('Time (s)')
            ylabel('Frequency')
            title('ECG signal')
            xlim([0,5])
        
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