tiempo = 2; % Sample signal for ~2 seconds
tic    
Fs = X; % Change sampling rate
t1 = toc;
y  = wavrecord(Fs*tiempo, Fs);
t2 = toc;
soundsc(y,Fs)
t = linspace(t1,t2,Fs*tiempo); % Time vector
figure('Name',['Fs: ', num2str(Fs),' Hz'])
    plot(t,y*1000)
    grid on
    xlabel('Tiempo (s)');
    xlim([t1,t2])
    ylabel('Amplitud (mV)');
    title(['Time: ',num2str(t2 - t1),' s.'])
    pause(1);
    
save AX y t1 t2 t Fs tiempo

wavplay(y,Fs)