function spect_analysis(settings,seq_nr,q1, q2, q3, q4, t, fig_nr,save_on_off)

     %Perform spectral analysis on the trajectory data
     L = length(t);
     
     Fs = 7500;
     
     %Calculating spectral power density of body positions
     
     NFFT = 2^nextpow2(L); % Next power of 2 from length of y
     Y1 = fft(q1,NFFT)/L;
     Y2 = fft(q2,NFFT)/L;
     Y3 = fft(q3,NFFT)/L;
     Y4 = fft(q4,NFFT)/L;
     f = Fs/2*linspace(0,1,NFFT/2+1);

     %Plot single-sided amplitude spectrum.
     figure(fig_nr)
     title('Spectral analysis quaternion')
     subplot(4,1,1); plot(f,2*abs(Y1(1:NFFT/2+1))) 
     title('q1')
     xlabel('Frequency in Hz')
     ylabel('|Amplitude|')
     subplot(4,1,2); plot(f,2*abs(Y2(1:NFFT/2+1))) 
     title('q2')
     xlabel('Frequency in Hz')
     ylabel('|Amplitude|')
     subplot(4,1,3); plot(f,2*abs(Y3(1:NFFT/2+1))) 
     title('q3')
     xlabel('Frequency in Hz')
     ylabel('|Amplitude|')
     subplot(4,1,4); plot(f,2*abs(Y4(1:NFFT/2+1))) 
     title('q4')
     xlabel('Frequency in Hz')
     ylabel('|Amplitude|')
     
     % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/spectral1'], 'fig')
    
    saveas(fig_nr, [char(settings.plot_folders(2)) '/spectral/spectral1_' int2str(seq_nr)], 'fig')
    
    end


end

