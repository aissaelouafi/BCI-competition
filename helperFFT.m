function helperFFT(bin, yVal,titleStr)
%Copyright 2014 The MathWorks, Inc
close all;clc;
figure;
plot(bin, yVal,'Color',[0,0,1],'LineWidth',1.5); box on; grid on;
xlabel('Frequency (Hz).'); 
if strcmp(titleStr,'Phase Response');
ylabel('Radians');
title('FFT - Phase Response');
else
ylabel('Magnitude');
title('FFT - Magnitude Response');
end