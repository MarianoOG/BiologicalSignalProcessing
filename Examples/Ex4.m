% Algorithm for a digital filter from the 60Hz cut off frequency low pass
% filter.

% H(z) = (0.3815*z + 0.3815)/(2.3815*z + 1.685);

% Signal:
fs = 1000;
t = 0:1/fs:1-1/fs;
fn = [10,30,60,80,100];
x = 0*t;
for i = 1:length(fn), x = x + sin(2*pi*fn(i)*t); end

% Filtered signal:
y = 0*t;
y(1) = x(1)*0.3815/2.3815;
for n = 2:length(x)
    y(n) = (y(n-1)*1.6185 + x(n)*0.3815 + x(n-1)*0.3815)/2.3815;
end

% FFT:
X = fft(x);
Y = fft(y);
Ie = floor(length(Y)/2);
f = (0:Ie-1)*((fs/2)/(Ie-1));

figure
    subplot(2,2,1)
        plot(t,x)
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Original time domain')
    subplot(2,2,3)
        plot(f,abs(X(1:Ie))/fs);
        title('Original frequency domain')
        xlabel('Frecuency (Hz)')
        ylabel('|X(f)|')
    subplot(2,2,2)
        plot(t,y)
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Filtered time domain')
    subplot(2,2,4)
        plot(f,abs(Y(1:Ie))/fs);
        title('Filtered frequency domain')
        xlabel('Frecuency (Hz)')
        ylabel('|X(f)|')