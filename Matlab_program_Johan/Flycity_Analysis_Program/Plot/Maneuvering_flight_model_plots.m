function Maneuvering_flight_model_plots( settings, pathDB )

    maneuver = load('FM_man_wb.mat');

    % plot a_dev vs Fx_forward:
    
    F_star_ax = maneuver.man_ax.F_star;
    
    N_ax = size(F_star_ax,2);
    
    a_dev_fx_eta_L = maneuver.man_ax.a_dev_eta_L;
    a_dev_fx_eta_R = maneuver.man_ax.a_dev_eta_R;
    
    a_dev_fx_eta = (a_dev_fx_eta_L+a_dev_fx_eta_R)./2;
    
    [Fx_star_ax, Fx_sort_id] = sort(F_star_ax(1,:));
    
    a_dev_fx_eta_sort = a_dev_fx_eta(:,Fx_sort_id);
    
    b_eta_n = pathDB.maneuver.c_fit_ax.b_eta_n(:,1,1);
    
    b_eta_p = pathDB.maneuver.c_fit_ax.b_eta_p(:,1,1);
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(2,3,1); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(1,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(1) 0],[Fx_star_ax(1)*b_eta_n(1) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(1)],'r')
    hold off
    ylabel('a_{\eta}(1)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,2); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(2,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(2) 0],[Fx_star_ax(2)*b_eta_n(2) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(2)],'r')
    hold off
    ylabel('a_{\eta}(2)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,3); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(3,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(3) 0],[Fx_star_ax(3)*b_eta_n(3) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(3)],'r')
    hold off
    ylabel('a_{\eta}(3)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,4); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(4,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(4) 0],[Fx_star_ax(4)*b_eta_n(4) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(4)],'r')
    hold off
    xlabel('F^*_{a x}')
    ylabel('a_{\eta}(4)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,5); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(5,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(5) 0],[Fx_star_ax(5)*b_eta_n(5) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(5)],'r')
    hold off
    xlabel('F^*_{a x}')
    ylabel('a_{\eta}(5)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,6); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(6,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(6) 0],[Fx_star_ax(6)*b_eta_n(6) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(6)],'r')
    hold off
    xlabel('F^*_{a x}')
    ylabel('a_{\eta}(6)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    hold off
    [~,h1] = suplabel('Deviation coefficients a_{\eta} downstroke vs F^*_{a x}', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Maneuvering_flight_model_plots/F_ax_down_vs_a_eta')],'fig')
   
    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(2,3,1); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(16,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(16) 0],[Fx_star_ax(16)*b_eta_n(16) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(16)],'r')
    hold off
    ylabel('a_{\eta}(16)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,2); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(17,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(17) 0],[Fx_star_ax(17)*b_eta_n(17) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(17)],'r')
    hold off
    ylabel('a_{\eta}(17)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,3); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(18,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(18) 0],[Fx_star_ax(18)*b_eta_n(18) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(18)],'r')
    hold off
    ylabel('a_{\eta}(18)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,4); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(19,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(19) 0],[Fx_star_ax(19)*b_eta_n(19) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(19)],'r')
    hold off
    xlabel('F^*_{a x}')
    ylabel('a_{\eta}(19)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,5); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(20,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(20) 0],[Fx_star_ax(20)*b_eta_n(20) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(20)],'r')
    hold off
    xlabel('F^*_{a x}')
    ylabel('a_{\eta}(20)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    subplot(2,3,6); hold on
    plot(Fx_star_ax,a_dev_fx_eta_sort(21,:),'o','MarkerSize',1.5,'MarkerEdgeColor', [0.5 0.5 0.5])
    plot([Fx_star_ax(21) 0],[Fx_star_ax(21)*b_eta_n(21) 0],'b')
    plot([0 Fx_star_ax(N_ax)],[0 Fx_star_ax(N_ax)*b_eta_p(21)],'r')
    hold off
    xlabel('F^*_{a x}')
    ylabel('a_{\eta}(21)')
    xlim([-4e-5 5e-5])
    ylim([-0.5 0.5])
    hold off
    [~,h1] = suplabel('Deviation coefficients a_{\eta} upstroke vs F^*_{a x}', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Maneuvering_flight_model_plots/F_ax_up_vs_a_eta')],'fig')
    
    % plot a_dev vs Mx:
    
    
    
    M_star_wx = maneuver.man_wx.M_star;

    N_wx = size(M_star_wx,2);
    
    a_dev_mx_theta_L = maneuver.man_wx.a_dev_theta_L;
    a_dev_mx_theta_R = maneuver.man_wx.a_dev_theta_R;
    
    Mx_star_wx = zeros(N_wx,1);
    
    a_dev_wx_theta_L = zeros(size(a_dev_mx_theta_L));
    a_dev_wx_theta_R = zeros(size(a_dev_mx_theta_R));
    
    for i = 1:N_wx
        
        if M_star_wx(1,i) < 0
        
            a_dev_wx_theta_L(:,i) = a_dev_mx_theta_R(:,i);
            a_dev_wx_theta_R(:,i) = a_dev_mx_theta_L(:,i);
            
            Mx_star_wx(i) = -M_star_wx(1,i);
        
        else
        
            a_dev_wx_theta_L(:,i) = a_dev_mx_theta_L(:,i);
            a_dev_wx_theta_R(:,i) = a_dev_mx_theta_R(:,i);
            
            Mx_star_wx(i) = M_star_wx(1,i);
            
        end
            
    end
    
    b_theta_L = pathDB.maneuver.c_fit_wx.b_theta_L(:,1,4);
    b_theta_R = pathDB.maneuver.c_fit_wx.b_theta_R(:,1,4);
    
    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(2,3,1); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(1,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(1,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(1)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(1)],'b')
    hold off
    ylabel('a_{\theta}(1)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,2); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(2,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(2,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(2)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(2)],'b')
    hold off
    ylabel('a_{\theta}(2)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,3); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(3,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(3,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(3)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(3)],'b')
    hold off
    ylabel('a_{\theta}(3)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,4); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(4,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(4,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(4)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(4)],'b')
    hold off
    xlabel('M^*_{\omega x}')
    ylabel('a_{\theta}(4)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,5); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(5,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(5,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(5)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(5)],'b')
    hold off
    xlabel('M^*_{\omega x}')
    ylabel('a_{\theta}(5)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,6); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(6,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(6,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(6)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(6)],'b')
    hold off
    xlabel('M^*_{\omega x}')
    ylabel('a_{\theta}(6)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    hold off
    [~,h1] = suplabel('Deviation coefficients a_{\theta L} and a_{\theta R} downstroke vs M^*_{\omega x}', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Maneuvering_flight_model_plots/M_wx_down_vs_a_theta')],'fig')
    
    figure(4)
    hFig = figure(4);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1200 800]);
    hold on
    subplot(2,3,1); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(14,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(14,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(14)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(14)],'b')
    hold off
    ylabel('a_{\theta}(14)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,2); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(15,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(15,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(15)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(15)],'b')
    hold off
    ylabel('a_{\theta}(15)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,3); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(16,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(16,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(16)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(16)],'b')
    hold off
    ylabel('a_{\theta}(16)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,4); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(17,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(17,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(17)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(17)],'b')
    hold off
    xlabel('M^*_{\omega x}')
    ylabel('a_{\theta}(17)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,5); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(18,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(18,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(18)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(18)],'b')
    hold off
    xlabel('M^*_{\omega x}')
    ylabel('a_{\theta}(18)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    subplot(2,3,6); hold on
    plot(Mx_star_wx,a_dev_wx_theta_L(19,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'r')
    plot(Mx_star_wx,a_dev_wx_theta_R(19,:),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_L(19)],'r')
    plot([0 max(Mx_star_wx)],[0 max(Mx_star_wx)*b_theta_R(19)],'b')
    hold off
    xlabel('M^*_{\omega x}')
    ylabel('a_{\theta}(19)')
    xlim([0 2.5e-5])
    ylim([-0.4 0.4])
    hold off
    [~,h1] = suplabel('Deviation coefficients a_{\theta L} and a_{\theta R} upstroke vs M^*_{\omega x}', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) char('/Maneuvering_flight_model_plots/M_wx_up_vs_a_theta')],'fig')
    

end

