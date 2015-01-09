function plot_alfa_beta(settings, seq_nr,alfa, beta, u_f, v_f, w_f, qb1_f, qb2_f, qb3_f, qb4_f, t, fig_nr1, save_on_off)
% 
%     U_body = zeros(3,length(u_f));
% 
%     for i = 1:length(u_f)
%         
%         U_vect = [u_f(i); v_f(i); w_f(i)];
%         
%         DCM = quat2matNEW([qb1_f(i) qb2_f(i) qb3_f(i) qb4_f(i)]);
%         
%         U_body(:,i) = DCM'*U_vect;
%         
%     end

    

    
    % Plot the angles alfa and beta w.r.t. the body.
    
    figure(fig_nr1)
    subplot(2,1,1); plot(t,radtodeg(alfa))
    title('Plot aoa and sideslip w.r.t. body.')
    xlabel('t [s]')
    ylabel('alfa [deg]')
    subplot(2,1,2); plot(t,radtodeg(beta))
    xlabel('t [s]')
    ylabel('beta [deg]')
    
%     figure(fig_nr2)
%     subplot(3,1,1); plot(t,U_body(1,:))
%     title('Body velocity components')
%     xlabel('t [s]')
%     ylabel('u [m/s]')
%     subplot(3,1,2); plot(t,U_body(2,:))
%     xlabel('t [s]')
%     ylabel('v [m/s]')
%     subplot(3,1,3); plot(t,U_body(3,:))
%     xlabel('t [s]')
%     ylabel('w [m/s]')
    
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/alfa_beta_body'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/alfa_beta_body/alfa_beta' int2str(seq_nr)], 'fig')
    
    end
    


end

