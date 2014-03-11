function [val_notch] = filter_borf_notch(val,cutoff,Q)

y = [val;val;val];
sampling = length(val); % sampling frequency = number of datapoints Hz  
R = cutoff/sampling;

%% butterworth
% % nyquist = sampling/2; % nyquist frequency is half of the sampling frequency  
% % Wn = cutoff/nyquist; % normalized cutoff frequency  
% % [b,a]=butter(n,Wn); % get Butterworth coefficients  
% [b,a]=butter(n,2*R,'low'); % get Butterworth coefficients, 2* because of nyquist
% y_butter=filtfilt(b,a,y); % apply filter to data  
% val_butter = y_butter(length(val)+1:2*length(val));
% 
% % plot
% figure
% freqz(b,a,128,sampling), title('filter characteristics')
% figure
% plot(val)
% hold on
% plot(val_butter,'r')

%% Notch
[b,a] = iirnotch(2*R,2*R/Q);
y_notch=filtfilt(b,a,y); % apply filter to data  
val_notch = y_notch(length(val)+1:2*length(val));
% fvtool(b,a);

% plot
figure
freqz(b,a,128,sampling), title('filter characteristics')
figure
plot(val)
hold on
plot(val_notch,'r')

%% lowpass
% b = R;
% a = [1 a1-1];
% y_lowpass = filtfilt(b,a,y);

