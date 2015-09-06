ecg_data = rdsamp('mitdb/100', 'begin', '17:20', 'maxt', '5:00');
t = ecg_data(:,1);
x = ecg_data(:,2);

Fs = 360;

%% Filtracja
y = sgolayfilt(x,3,1001);
c = smooth(x,500,'moving');

subplot(3,1,1), plot(t,x), hold on, plot(t,y,'r','linewidth',2), plot(t,c,'g','linewidth',2),
title('ECG + wandering signal'); 
legend('ECG', 'Savitzky-Golay', 'moving average');

subplot(3,1,2), plot(t,x-y),
title('ECG with wandering filtered out using Savitzky-Golay filter');
subplot(3,1,3), plot(t,x-c),
title('ECG with wandering filtered out using moving average');
