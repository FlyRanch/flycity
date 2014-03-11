function [y_butter] = filter_borf_raw(y,cutoff,n,n_wb)

% cutoff = 20; % cutoff frequency = 20 Hz  
% n = 3; %filter order (higher is sharper cutoff)  

sampling = length(y)/n_wb; % sampling frequency = number of datapoints Hz  
R = cutoff/sampling;

%% butterworth
% nyquist = sampling/2; % nyquist frequency is half of the sampling frequency  
% Wn = cutoff/nyquist; % normalized cutoff frequency  
% [b,a]=butter(n,Wn); % get Butterworth coefficients  
[b,a]=butter(n,2*R,'low'); % get Butterworth coefficients, 2* because of nyquist
y_butter=filtfilt(b,a,y); % apply filter to data  


%% plot
% figure
% freqz(b,a,128,sampling), title('low pass filter characteristics')
% figure
% plot(y)
% hold on
% plot(y_butter,'r')

