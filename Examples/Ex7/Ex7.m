function [] = Ex7()
    
%     H = [9,15,25,14,10,18,0,16,5,19,16,20];
%     M = [39,56,93,61,50,75,32,85,42,70,66,80];
%     C = cov(H,M)
%     [A,B] = eig(C)
%     [A,B] = pcacov([H',M']);
    [time,signal,Fs,~] = readPhysionet('ecg1');
    n = length(signal)/10;
    A = zeros(n,10);
    for i = 1:10
        A(:,i) = signal((i-1)*n+1:i*n);
    end
    pcacov(A)
    figure, plot(A)     
    figure
            plot(time,signal)    
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('ECG signal')
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