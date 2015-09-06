butterworth = importdata('butterworth.txt')';
savitzkygolay = importdata('savitzkygolay.txt')';
signal = importdata('signal_base.txt')';
wavelet = importdata('wavelet.txt')';
t = 1:1/360:14400/360+1;
t = t(1:end-1);
plot(t,butterworth,t,wavelet,t,savitzkygolay)
legend('butterworth','wavelet','savitzkygolay');