function [val_butter] = filter_borf_butter_dt(val,cutoff,n,dt)

% cutoff = 20; % cutoff frequency = 20 Hz  
% n = 3; %filter order (higher is sharper cutoff)  

y = [val;val;val];
sampling = length(val); % sampling frequency = number of datapoints Hz  
R = cutoff/sampling;

%% butterworth
% nyquist = sampling/2; % nyquist frequency is half of the sampling frequency  
% Wn = cutoff/nyquist; % normalized cutoff frequency  
% [b,a]=butter(n,Wn); % get Butterworth coefficients  
[b,a]=butter(n,2*R,'low'); % get Butterworth coefficients, 2* because of nyquist
y_butter=filtfilt(b,a,y); % apply filter to data  

%% lowpass
% b = R;
% a = [1 a1-1];
% y_lowpass = filtfilt(b,a,y);



val_butter = y_butter(length(val)+1:2*length(val));


%% plot
% figure
% freqz(b,a,128,sampling), title('low pass filter characteristics')
% figure
% plot(val)
% hold on
% plot(val_butter,'r')

