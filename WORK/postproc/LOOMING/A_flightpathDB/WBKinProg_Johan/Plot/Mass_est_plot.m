function Mass_est_plot(settings,pathDB,seq_nr,save_on_off, fig_nr)


    % Program that plots a bar graph from the estimated mass based on
    % measured acceleration and obtained force by the quasi-steady state
    % model.

    
    for i = 1:find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0,1, 'last')
        
        a = 1:find(isnan(pathDB.wingbeat_time(i,:,seq_nr))==0,1, 'last');
        
        b = pathDB.wingbeat_time(i,a,seq_nr);
        
        m_bar(i) = mean(pathDB.m_est_point(b,seq_nr));
                
    end
    
    
    figure(fig_nr)
    bar(m_bar)
    title('Mass estimation by quasi-steady model')
    xlabel('wingbeat nr')
    ylabel('M_{est}/M_{winglength}')
    
        % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/M_est'], 'fig')
    
    end

end
