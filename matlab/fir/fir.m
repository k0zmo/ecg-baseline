%% Pobranie danych
ecg_data = rdsamp('mitdb/101', 'maxt', ':40');
t = ecg_data(:,1);
x = ecg_data(:,2);

% Hhp: filterbuilder, IIR, 10-th order, 3dB point, Fs = 360, F3dB = 0.5Hz,
% Butterworth, Direct-form I SOS
yzp = sosfiltfilt(Hhp.sosMatrix, x) * Hhp.ScaleValues(end);
%y = sosfilt(Hhp.sosMatrix, x) * Hhp.ScaleValues(end);
y = sosfilt1(Hhp.sosMatrix, x) * Hhp.ScaleValues(end);

plot(linspace(0,1,length(x)),x,'r'), hold on
plot(linspace(0,1,length(yzp)),yzp,'linewidth',1)
plot(linspace(0,1,length(x)),y,'k');
legend('sygnal oryginalny','zero-phase sosfilt','standard sosfilt');