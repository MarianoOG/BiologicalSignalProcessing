Fs = 1000;
tfin = 3;
t = 0:1/Fs:tfin;
t2 = -tfin:1/Fs:tfin;
x = 0*t;
x(1:Fs*2+1) = ones(1,Fs*2+1);
y = 0*t;
y(1:Fs+1) = 1:-1/Fs:0;

Rxy = xcorr(x,y)/Fs;
Ryx = xcorr(y,x)/Fs;

figure, plot(t,x,t,y,t2,Rxy,t2,Ryx)
	legend('x(t)','y(t)','Rxy(t)','Ryx(t)',2)
    xlabel('Time (s)')
    ylabel('Amplitude')
    title('Correlation')
    grid on

R = 0*t2;
R(Fs*2+1:Fs*3+1) = (t2(Fs*2+1:Fs*3+1).^2+2*t2(Fs*2+1:Fs*3+1)+1)/2;
R(Fs*3+1:Fs*4+1) = 0.5*ones(1,Fs+1);
R(Fs*4+1:Fs*5+1) = -t2(Fs*4+1:Fs*5+1).^2./2+t2(Fs*4+1:Fs*5+1);

figure
    plot(t2,R)
    xlabel('Time (s)')
    ylabel('Amplitude')
    title('Correlation')
    grid on
    