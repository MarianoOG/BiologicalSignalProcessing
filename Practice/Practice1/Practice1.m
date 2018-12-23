%% Practice 1 -  Sampling theorem
%% 1. Generate a sine wave at 50Hz using the following sampling frequencies:
VFs = [1000,900,800,700,600,500,400,300,200,100,50];
t = {};
x = {};
for i = 1:length(VFs);
    Fs = VFs(i);
    t{i} = 0:1/Fs:0.08;
    x{i} = sin(2*pi*50*t{i});
    if 1 == mod(i,4), figure, end
    if 0 == mod(i,4)
        subplot(4,1,4)
    else
        subplot(4,1,mod(i,4))
    end
        plot(t{i},x{i})
        title(['50 Hz sine wave at fs = ',num2str(Fs)])
        xlabel('Time (s)')
        ylabel('A(v)')
end
%% 2 Plot a segment of your voice recordings for each sampling frequency:
Rec = {};
load A5000
Rec{1,1} = y;
Rec{1,2} = t;
Rec{1,3} = Fs;
load A11000
Rec{2,1} = y;
Rec{2,2} = t;
Rec{2,3} = Fs;
load A22000
Rec{3,1} = y;
Rec{3,2} = t;
Rec{3,3} = Fs;
load A44100
Rec{4,1} = y;
Rec{4,2} = t;
Rec{4,3} = Fs;
for i = 1:4
    y = Rec{i,1};
    t = Rec{i,2};
    Fs = Rec{i,3};
    y = y(0.7*Fs:1.7*Fs);
    t = t(0.7*Fs:1.7*Fs);
    Y = fft(y)/Fs;
    Ie = floor(length(Y)/2);
    f = (0:Ie-1)*((Fs/2)/(Ie-1));
    % Graficar:
    if mod(i,2) == 1, figure, end
    j = i;
    if i > 2; j = j - 2; end
    subplot(2,2,j*2-1)
        plot(t(0.85*Fs:0.86*Fs),y(0.85*Fs:0.86*Fs))
        xlabel('Tiempo (s)');
        xlim([t(0.85*Fs),t(0.86*Fs)])
        ylabel('Amplitud (mV)');
        title(['Voz a fs = ',num2str(Fs)])
    subplot(2,2,j*2)
        plot(f,abs(Y(1:Ie))/length(Y));
        title(['Frecuency domain at Fs = ',num2str(Fs)])
        xlabel('Frecuency (Hz)')
        xlim([0,2000])
        ylabel('|X(f)|')
end