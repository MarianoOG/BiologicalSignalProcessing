function [] = Practice2(  )
% Calculate the FT of the following functions:
    Data = {};
% 1) x(t) = A*sin(2*pi*f0*t); f0 = 50 Hz; Fs = 1000 Hz;
    Fs = 1000;
    f0 = 50;
    t = 0:1/Fs:0.2;
    x = sin(2*pi*f0*t);
    [f,X] = fourierVector(x,Fs);
    Data = save(t,x,f,X,Fs,Data);
% 2) x(t) = A*rect(t/T); Fs = 100 Hz;
    Fs = 100;
    T = 2;
    [t,x] = rect(T,Fs);
    [f,X] = fourierVector(x,Fs);
    Data = save(t,x,f,X,Fs,Data);
% 3) x(t) = A*exp(-a*t); Fs = 10 Hz;
    Fs = 10;
    a = 5;
    t = 0:1/Fs:4;
    x = exp(-a*t);
    [f,X] = fourierVector(x,Fs);
    Data = save(t,x,f,X,Fs,Data);
% 4) x(t) = A*sin() + sin(); f01 = 80 Hz; f02 = 120 Hz; Fs = propose;
    Fs = 1200;
    t = 0:1/Fs:0.1;
    f1 = 80;
    f2 = 120;
    x = 2*sin(2*pi*f1*t) + sin(2*pi*f2*t);
    [f,X] = fourierVector(x,Fs);
    Data = save(t,x,f,X,Fs,Data);
% 5) x(t) = A*tri(t/T); Fs = 100 Hz;
    Fs = 100;
    T = 1;
    [t,x] = tri(T,Fs);
    [f,X] = fourierVector(x,Fs);
    Data = save(t,x,f,X,Fs,Data);
% Graphics:
    limx = [150,20,5,400,5];
    for i = 1:size(Data,1)
        figure
            subplot(2,1,1)
                plot(Data{i,1},Data{i,2})
                xlabel('Time (s)')
                ylabel('Amplitude (v)')
                ylim([1.1*min(Data{i,2}),1.1*max(Data{i,2})])
                title(['Time domain at Fs = ',num2str(Data{i,5})])
            subplot(2,1,2)
                plot(Data{i,3},Data{i,4})
                xlabel('Frequency (Hz)')
                xlim([0,limx(i)])
                ylabel('Amplitude (v)')
                ylim([0,1.1*max(Data{i,4})])
                title(['Frequency domain at Fs = ',num2str(Data{i,5})])
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

function [Data] = save(t,x,f,X,Fs,Data)
% Stores data in a cell
    l = size(Data,1);
    Data{l+1,1} = t;
    Data{l+1,2} = x;
    Data{l+1,3} = f;
    Data{l+1,4} = X;
    Data{l+1,5} = Fs;
end

function [t,x] = rect(T,Fs)
% Creates a rect function
    t = -2.5*T:1/Fs:2.5*T;
    x = 0*t;
    cero = floor(length(x)/2)+1;
    x(cero-floor(T/2)*Fs+1:cero+floor(T/2)*Fs) = ones(1,T*Fs);
end

function [t,x] = tri(T,Fs)
% Creates a rect function
    T = 2*T;
    t = -2.5*T:1/Fs:2.5*T;
    x = 0*t;
    cero = floor(length(x)/2)+1;
    x(cero-floor(T/2)*Fs+1:cero) = (0:Fs*floor(T/2)-1)/(floor(T/2)*Fs-1);
    x(cero:cero+floor(T/2)*Fs) = 1-(0:floor(T/2)*Fs)/(floor(T/2)*Fs);
end