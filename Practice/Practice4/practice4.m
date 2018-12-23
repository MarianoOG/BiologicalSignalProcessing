function [] = practice4 ()
% Practice 4: Convolution and correlation

%% Part I: convolution
%   1) the 3 exersices seen in class
%   2) Sketch the original signals in MATLAB
%   3) Solve the convolution using the 'conv' command
%   4) Solve the convolution using the 
%   5) Sketch the solution in MATLAB to compare with your solution.
	
    % Variables:
    Fs = 1000;
    tfin = 10;
    t = 0:1/Fs:tfin;
    t2 = 0:1/Fs:2*tfin;
    funciones = cell(3,4);
    
    % 1)
    x = 0*t;
    x(4*Fs+1:6*Fs+1) = 2*ones(1,Fs*2+1);
    h = 0*t;
    h(1:4*Fs+1) = ones(1,4*Fs+1);
    y = conv(x,h)/Fs;
    funciones = salvar(x,h,y,funciones,1);
    
    % 2)
    x = exp(-t);
    h = 0*t;
    h(1:Fs+1) = 0:1/Fs:1;
    y = conv(x,h)/Fs;
    funciones = salvar(x,h,y,funciones,2);
    
    % 3)
    x = 2*exp(-t);
    h = exp(-2*t);
    y = conv(x,h)/Fs;
    funciones = salvar(x,h,y,funciones,3);
    
    % 4)
    x = 0*t;
    x(2*Fs+1:3*Fs+1) = 2*(0:1/Fs:1);
    x(3*Fs+1:4*Fs+1) = 2*(1:-1/Fs:0);
    h = exp(-t);
    y = conv(x,h)/Fs;
    funciones = salvar(x,h,y,funciones,4);
        
    % Graficar:
    lim = [11,6,6,9];
    for i = 1:size(funciones,2)
        figure, plot(t,funciones{1,i},t,funciones{2,i},t2,funciones{3,i})
            legend('x(t)','h(t)','y(t)')
            xlabel('Time (s)')
            ylabel('Amplitude')
            title(['Convolution - Exercise ',num2str(i)])
            xlim([0,lim(i)])
    end
    
    clear all
%% Part II: Correlation
% Given x(t) and y(t) find: a) Rxx, b) Ryy, c) Rxy, d) Ryx

% 1) Find the analytic solution.
% 2) Implement the correlation using the command "xcorr"
% 3) Compare the graphcs of your analytic and the "xcorr" solution
% 4) Select a biosignal and use correlation to obtain information comment
%    your results

    Fs = 1000;
    tlim = 2;
    t = -tlim:1/Fs:tlim;
    t2 = -2*tlim:1/Fs:2*tlim;
    funciones = cell(6,1);
    
    % Ex 1)
    x = 0*t;
    x(1*Fs+1:2*Fs+1) = ones(1,Fs+1);
    y = 0*t;
    y(2*Fs+1:3*Fs+1) = t(2*Fs+1:3*Fs+1);
    funciones{1,1} = x;
    funciones{2,1} = y;
    funciones{3,1} = xcorr(x,x)/Fs; % Rxx
    funciones{4,1} = xcorr(x,y)/Fs; % Rxy
    funciones{5,1} = xcorr(y,x)/Fs; % Ryx
    funciones{6,1} = xcorr(y,y)/Fs; % Ryy
    
    % Graficar:
    lim = [2.5];
    for i = 1:size(funciones,2)
        figure, plot(t,funciones{1,i},t,funciones{2,i},t2,funciones{3,i},...
            t2,funciones{4,i},t2,funciones{5,i},t2,funciones{6,i})
            legend('x(t)','y(t)','Rxx','Rxy','Ryx','Ryy')
            xlabel('Time (s)')
            ylabel('Amplitude')
            title(['Correlation - Exercise ',num2str(i)])
            xlim([-lim(i),lim(i)])
    end
    
    % Biosignal:
    [time,signal,Fs,~] = readPhysionet('ecg1');
    QRS = 0*time;
    QRS(11:18) = signal(68:75); % Extrae conjunto un conjunto QRS
    time2 = -max(time)+1/Fs:1/Fs:max(time)-1/Fs;
    R = xcorr(signal,QRS)/Fs;
    Ra = (R(length(time):end)>0).*R(length(time):end)/max(R);
    figure, 
    subplot(2,2,1)
        plot(time,signal)
        xlabel('Time (s)')
        ylabel('Amplitude')
        title('ECG signal')
        grid on
        xlim([0,3])
    subplot(2,2,2)
        plot(time,QRS)
        xlabel('Time (s)')
        ylabel('Amplitude')
        title('Extracted QRS from signal')
        grid on
        xlim([0,0.3])
    subplot(2,1,2)
        plot(time2,R)
        xlabel('Time (s)')
        ylabel('Amplitude')
        title('Correlation')
        grid on
        xlim([0,3])
    
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

function [c] = salvar(x,h,y,c,n)
% Save vector:
    c{1,n} = x;
    c{2,n} = h;
    c{3,n} = y;
end