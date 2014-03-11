function [low_data] = low_pass_filter(data,f_cutoff,n)

%% Implemention of a Low pass Butterworth filter using the filtfilt command

%% change the data, sampling frequency, cutoff frequency and order of the filter based on your requirements
%% I have provided sample data file- a noisy sine wave

% File written by Avinash Parnandi. http://robotics.usc.edu/~parnandi

% xx = [1:100];
% data = sin(188*xx)+rand(1,100);  % noisy data; change it to whatever is your data
% 
% f=30;%  sampling frequency
% f_cutoff = 5; % cutoff frequency

f=length(data);%  sampling frequency
fnorm =f_cutoff/(f/2); % normalized cut off freq, you can change it to any value depending on your requirements

[b1,a1] = butter(n,fnorm,'low'); % Low pass Butterworth filter of order n
low_data = filtfilt(b1,a1,data); % filtering

figure
freqz(b1,a1,128,f), title('low pass filter characteristics')
figure
% subplot(2,1,1), plot(data), title('Actual data')
% subplot(2,1,2), plot(low_data), title('Filtered data')
plot(data)
hold on
plot(low_data,'r')

