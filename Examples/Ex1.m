close all; clc; clear all;
%% Excersice 1:
% Generate a sine wave
Fs = 1000;
f0 = 50;
t = 0:1/Fs:0.2;
x = sin(2*pi*f0*t);
X = fft(x);
Ie = floor(length(X)/2);
f = (0:Ie-1)*((Fs/2)/(Ie-1));
% Plot cosine wave, plot FT of cosine
figure
    subplot(2,1,1)
        plot(t,x)
        title('Time domain')
        xlabel('Time (s)')
        ylabel('x(t)')
    subplot(2,1,2)
        stem(f,abs(X(1:Ie))/length(X));
        title('Frecuency domain')
        xlabel('Frecuency (Hz)')
        ylabel('|X(f)|')
% Generate a cosine wave
Fs = 1000;
f0 = 50;
t = 0:1/Fs:0.2;
x = cos(2*pi*f0*t);
X = fft(x);
Ie = floor(length(X)/2);
f = (0:Ie-1)*((Fs/2)/(Ie-1));
% Plot cosine wave, plot FT of cosine
figure
    subplot(2,1,1)
        plot(t,x)
        title('Time domain')
        xlabel('Time (s)')
        ylabel('x(t)')
    subplot(2,1,2)
        stem(f,abs(X(1:Ie))/length(X));
        title('Frecuency domain')
        xlabel('Frecuency (Hz)')
        ylabel('|X(f)|')
