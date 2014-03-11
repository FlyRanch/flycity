function [n] = plot_dyn_sim( sim_data, case_nr, n )

    N               = length(sim_data.T);

    T               = sim_data.T;
    xyz             = sim_data.xyz;
    qb              = sim_data.qb;
    vb              = sim_data.vb;
    wb              = sim_data.wb;
    ab              = sim_data.ab;
    w_dot_b         = sim_data.w_dot_b;
    vb_mean         = sim_data.vb_mean;
    wb_mean         = sim_data.wb_mean;
    ab_mean         = sim_data.ab_mean;
    w_dot_b_mean    = sim_data.w_dot_b_mean;
    FI_vel_b        = sim_data.FI_vel_b;
    MI_vel_b        = sim_data.MI_vel_b;
    FI_vel_strk     = sim_data.FI_vel_strk;
    MI_vel_strk     = sim_data.MI_vel_strk;
    FI_acc_b        = sim_data.FI_acc_b;
    MI_acc_b        = sim_data.MI_acc_b;
    FI_acc_strk     = sim_data.FI_acc_strk;
    MI_acc_strk     = sim_data.MI_acc_strk;
    LM              = sim_data.LM;
    AM              = sim_data.AM;
    FI_vel_b_mean   = sim_data.FI_vel_b_mean;
    MI_vel_b_mean   = sim_data.MI_vel_b_mean;
    FI_vel_strk_mean= sim_data.FI_vel_strk_mean;
    MI_vel_strk_mean= sim_data.MI_vel_strk_mean;
    FI_acc_b_mean   = sim_data.FI_acc_b_mean;
    MI_acc_b_mean   = sim_data.MI_acc_b_mean;
    FI_acc_strk_mean= sim_data.FI_acc_strk_mean;
    MI_acc_strk_mean= sim_data.MI_acc_strk_mean;
    LM_mean         = sim_data.LM_mean;
    AM_mean         = sim_data.AM_mean;
    KE              = sim_data.KE;
    KE_lin          = sim_data.KE_lin;
    KE_ang          = sim_data.KE_ang;
    FA_b            = sim_data.FA_b;
    MA_b            = sim_data.MA_b;
    FA_strk         = sim_data.FA_strk;
    MA_strk         = sim_data.MA_strk;
    alfa_L          = sim_data.alfa_L;
    alfa_R          = sim_data.alfa_R;
    alfa_dot_L      = sim_data.alfa_dot_L;
    alfa_dot_R      = sim_data.alfa_dot_R;
    FA_b_mean       = sim_data.FA_b_mean;
    MA_b_mean       = sim_data.MA_b_mean;
    FA_strk_mean    = sim_data.FA_strk_mean;
    MA_strk_mean    = sim_data.MA_strk_mean;
    Fg_b            = sim_data.Fg_b;
    Fg_strk         = sim_data.Fg_strk;
    
    figure(n)
    hFig = figure(n);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(3,1,1); plot(T,xyz(1,:))
    xlabel('t [s]')
    ylabel('x [mm]')
    subplot(3,1,2); plot(T,xyz(2,:))
    xlabel('t [s]')
    ylabel('y [mm]')
    subplot(3,1,3); plot(T,xyz(3,:))
    xlabel('t [s]')
    ylabel('z [mm]')
    hold off
    [~,h1] = suplabel('Position c.g. body', 't');
    set(h1,'FontSize',10)
    %print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FX_forward.eps'] ,'-depsc2');
    
    n = n+1;
    
    figure(n)
    hFig = figure(n);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(4,1,1); plot(T,qb(1,:))
    xlabel('t [s]')
    ylabel('qb_1')
    subplot(4,1,2); plot(T,qb(2,:))
    xlabel('t [s]')
    ylabel('qb_2')
    subplot(4,1,3); plot(T,qb(3,:))
    xlabel('t [s]')
    ylabel('qb_3')
    subplot(4,1,4); plot(T,qb(4,:))
    xlabel('t [s]')
    ylabel('qb_4')
    hold off
    [~,h1] = suplabel('Body quaternion', 't');
    set(h1,'FontSize',10)

    n = n+1;

    figure(n)
    hFig = figure(n);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(3,1,1); plot(T,vb(1,:),'r',T,ones(N,1)*vb_mean(1),'k')
    xlabel('t [s]')
    ylabel('v_x [mm/s]')
    subplot(3,1,2); plot(T,vb(2,:),'r',T,ones(N,1)*vb_mean(2),'k')
    xlabel('t [s]')
    ylabel('v_y [mm/s]')
    subplot(3,1,3); plot(T,vb(3,:),'r',T,ones(N,1)*vb_mean(3),'k')
    xlabel('t [s]')
    ylabel('v_z [mm/s]')
    hold off
    [~,h1] = suplabel('Velocity body (body ref. frame)', 't');
    set(h1,'FontSize',10)
    
    n = n+1;

    figure(n)
    hFig = figure(n);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(3,1,1); plot(T,wb(1,:),'r',T,ones(N,1)*wb_mean(1),'k')
    xlabel('t [s]')
    ylabel('\omega_x [rad/s]')
    subplot(3,1,2); plot(T,wb(2,:),'r',T,ones(N,1)*wb_mean(2),'k')
    xlabel('t [s]')
    ylabel('\omega_y [rad/s]')
    subplot(3,1,3); plot(T,wb(3,:),'r',T,ones(N,1)*wb_mean(3),'k')
    xlabel('t [s]')
    ylabel('\omega_z [rad/s]')
    hold off
    [~,h1] = suplabel('Angular velocity body (body ref. frame)', 't');
    set(h1,'FontSize',10)
    
    n = n+1;

    figure(n)
    hFig = figure(n);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(3,1,1); plot(T,ab(1,:),'r',T,ones(N,1)*ab_mean(1),'k')
    xlabel('t [s]')
    ylabel('a_x [mm/s^2]')
    subplot(3,1,2); plot(T,ab(2,:),'r',T,ones(N,1)*ab_mean(2),'k')
    xlabel('t [s]')
    ylabel('a_y [mm/s^2]')
    subplot(3,1,3); plot(T,ab(3,:),'r',T,ones(N,1)*ab_mean(3),'k')
    xlabel('t [s]')
    ylabel('a_z [mm/s^2]')
    hold off
    [~,h1] = suplabel('Acceleration body (body ref. frame)', 't');
    set(h1,'FontSize',10)

    n = n+1;

    figure(n)
    hFig = figure(n);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(3,1,1); plot(T,w_dot_b(1,:),'r',T,ones(N,1)*w_dot_b_mean(1),'k')
    xlabel('t [s]')
    ylabel('$\dot{\omega} x [rad/s^2]$','interpreter','latex')
    subplot(3,1,2); plot(T,w_dot_b(2,:),'r',T,ones(N,1)*w_dot_b_mean(2),'k')
    xlabel('t [s]')
    ylabel('$\dot{\omega} y [rad/s^2]$','interpreter','latex')
    subplot(3,1,3); plot(T,w_dot_b(3,:),'r',T,ones(N,1)*w_dot_b_mean(3),'k')
    xlabel('t [s]')
    ylabel('$\dot{\omega} z [rad/s^2]$','interpreter','latex')
    hold off
    [~,h1] = suplabel('Angular acceleration body (body ref. frame)', 't');
    set(h1,'FontSize',10)
    
    if case_nr == 0
    
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,FI_vel_b(1,:),'r',T,FI_acc_b(1,:),'b',T,ones(N,1)*FI_vel_b_mean(1),'r',T,ones(N,1)*FI_acc_b_mean(1),'b')
        xlabel('t [s]')
        ylabel('F_x')
        legend('F velocity','F acceleration','F vel mean','F acc mean')
        subplot(3,1,2); plot(T,FI_vel_b(2,:),'r',T,FI_acc_b(2,:),'b',T,ones(N,1)*FI_vel_b_mean(2),'r',T,ones(N,1)*FI_acc_b_mean(2),'b')
        xlabel('t [s]')
        ylabel('F_y')
        legend('F velocity','F acceleration','F vel mean','F acc mean')
        subplot(3,1,3); plot(T,FI_vel_b(3,:),'r',T,FI_acc_b(3,:),'b',T,ones(N,1)*FI_vel_b_mean(3),'r',T,ones(N,1)*FI_acc_b_mean(3),'b')
        xlabel('t [s]')
        ylabel('F_z')
        legend('F velocity','F acceleration','F vel mean','F acc mean')
        hold off
        [~,h1] = suplabel('Inertia forces (body ref. frame)', 't');
        set(h1,'FontSize',10)    

        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,MI_vel_b(1,:),'r',T,MI_acc_b(1,:),'b',T,ones(N,1)*MI_vel_b_mean(1),'r',T,ones(N,1)*MI_acc_b_mean(1),'b')
        xlabel('t [s]')
        ylabel('M_x')
        legend('M velocity','M acceleration')
        subplot(3,1,2); plot(T,MI_vel_b(2,:),'r',T,MI_acc_b(2,:),'b',T,ones(N,1)*MI_vel_b_mean(2),'r',T,ones(N,1)*MI_acc_b_mean(2),'b')
        xlabel('t [s]')
        ylabel('M_y')
        legend('M velocity','M acceleration')
        subplot(3,1,3); plot(T,MI_vel_b(3,:),'r',T,MI_acc_b(3,:),'b',T,ones(N,1)*MI_vel_b_mean(3),'r',T,ones(N,1)*MI_acc_b_mean(3),'b')
        xlabel('t [s]')
        ylabel('M_z')
        legend('M velocity','M acceleration')
        hold off
        [~,h1] = suplabel('Inertia moments (body ref. frame)', 't');
        set(h1,'FontSize',10)    

        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,FI_vel_strk(1,:),'r',T,FI_acc_strk(1,:),'b',T,ones(N,1)*FI_vel_strk_mean(1),'r',T,ones(N,1)*FI_acc_strk_mean(1),'b')
        xlabel('t [s]')
        ylabel('F_x')
        legend('F velocity','F acceleration')
        subplot(3,1,2); plot(T,FI_vel_strk(2,:),'r',T,FI_acc_strk(2,:),'b',T,ones(N,1)*FI_vel_strk_mean(2),'r',T,ones(N,1)*FI_acc_strk_mean(2),'b')
        xlabel('t [s]')
        ylabel('F_y')
        legend('F velocity','F acceleration')
        subplot(3,1,3); plot(T,FI_vel_strk(3,:),'r',T,FI_acc_strk(3,:),'b',T,ones(N,1)*FI_vel_strk_mean(3),'r',T,ones(N,1)*FI_acc_strk_mean(3),'b')
        xlabel('t [s]')
        ylabel('F_z')
        legend('F velocity','F acceleration')
        hold off
        [~,h1] = suplabel('Inertia forces (strokeplane ref. frame)', 't');
        set(h1,'FontSize',10)    

        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,MI_vel_strk(1,:),'r',T,MI_acc_strk(1,:),'b',T,ones(N,1)*MI_vel_strk_mean(1),'r',T,ones(N,1)*MI_acc_strk_mean(1),'b')
        xlabel('t [s]')
        ylabel('M_x')
        legend('M velocity','M acceleration')
        subplot(3,1,2); plot(T,MI_vel_strk(2,:),'r',T,MI_acc_strk(2,:),'b',T,ones(N,1)*MI_vel_strk_mean(2),'r',T,ones(N,1)*MI_acc_strk_mean(2),'b')
        xlabel('t [s]')
        ylabel('M_y')
        legend('M velocity','M acceleration')
        subplot(3,1,3); plot(T,MI_vel_strk(3,:),'r',T,MI_acc_strk(3,:),'b',T,ones(N,1)*MI_vel_strk_mean(3),'r',T,ones(N,1)*MI_acc_strk_mean(3),'b')
        xlabel('t [s]')
        ylabel('M_z')
        legend('M velocity','M acceleration')
        hold off
        [~,h1] = suplabel('Inertia moments (strokeplane ref. frame)', 't');
        set(h1,'FontSize',10)    

        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,LM(1,:),'r',T,ones(N,1)*LM_mean(1),'k')
        xlabel('t [s]')
        ylabel('LM_x')
        subplot(3,1,2); plot(T,LM(2,:),'r',T,ones(N,1)*LM_mean(2),'k')
        xlabel('t [s]')
        ylabel('LM_y')
        subplot(3,1,3); plot(T,LM(3,:),'r',T,ones(N,1)*LM_mean(3),'k')
        xlabel('t [s]')
        ylabel('LM_z')
        hold off
        [~,h1] = suplabel('Linear momentum (body ref. frame)', 't');
        set(h1,'FontSize',10)

        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,AM(1,:),'r',T,ones(N,1)*AM_mean(1),'k')
        xlabel('t [s]')
        ylabel('AM_x')
        subplot(3,1,2); plot(T,AM(2,:),'r',T,ones(N,1)*AM_mean(2),'k')
        xlabel('t [s]')
        ylabel('AM_y')
        subplot(3,1,3); plot(T,AM(3,:),'r',T,ones(N,1)*AM_mean(3),'k')
        xlabel('t [s]')
        ylabel('AM_z')
        hold off
        [~,h1] = suplabel('Angular momentum (body ref. frame)', 't');
        set(h1,'FontSize',10)

        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);  
        hold on
        plot(T,KE,'k')
        plot(T,KE_lin,'b')
        plot(T,KE_ang,'r')
        hold off
        title('Kinetic energy')
        xlabel('t [s]')
        ylabel('Energy [*1e-6 J]')
        legend('Total kin. E','Linear kin. E','Angular kin. E')
    
    else
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,FA_b(1,:),'b',T,ones(N,1)*FA_b_mean(1),'k')
        xlabel('t [s]')
        ylabel('F_x')
        subplot(3,1,2); plot(T,FA_b(2,:),'b',T,ones(N,1)*FA_b_mean(2),'k')
        xlabel('t [s]')
        ylabel('F_y')
        subplot(3,1,3); plot(T,FA_b(3,:),'b',T,ones(N,1)*FA_b_mean(3),'k')
        xlabel('t [s]')
        ylabel('F_z')
        hold off
        [~,h1] = suplabel('Aerodynamic forces (body ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,MA_b(1,:),'b',T,ones(N,1)*MA_b_mean(1),'k')
        xlabel('t [s]')
        ylabel('M_x')
        subplot(3,1,2); plot(T,MA_b(2,:),'b',T,ones(N,1)*MA_b_mean(2),'k')
        xlabel('t [s]')
        ylabel('M_y')
        subplot(3,1,3); plot(T,MA_b(3,:),'b',T,ones(N,1)*MA_b_mean(3),'k')
        xlabel('t [s]')
        ylabel('M_z')
        hold off
        [~,h1] = suplabel('Aerodynamic moments (body ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,FA_strk(1,:),'b',T,ones(N,1)*FA_strk_mean(1),'k')
        xlabel('t [s]')
        ylabel('F_x')
        subplot(3,1,2); plot(T,FA_strk(2,:),'b',T,ones(N,1)*FA_strk_mean(2),'k')
        xlabel('t [s]')
        ylabel('F_y')
        subplot(3,1,3); plot(T,FA_strk(3,:),'b',T,ones(N,1)*FA_strk_mean(3),'k')
        xlabel('t [s]')
        ylabel('F_z')
        hold off
        [~,h1] = suplabel('Aerodynamic forces (strokeplane ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,MA_strk(1,:),'b',T,ones(N,1)*MA_strk_mean(1),'k')
        ylabel('M_x')
        subplot(3,1,2); plot(T,MA_strk(2,:),'b',T,ones(N,1)*MA_strk_mean(2),'k')
        xlabel('t [s]')
        ylabel('M_y')
        subplot(3,1,3); plot(T,MA_strk(3,:),'b',T,ones(N,1)*MA_strk_mean(3),'k')
        xlabel('t [s]')
        ylabel('M_z')
        hold off
        [~,h1] = suplabel('Aerodynamic moments (strokeplane ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

%         figure(n)
%         hFig = figure(n);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 600 600]);
%         plot(T,radtodeg(alfa_L(end,:)),'r',T,radtodeg(alfa_R(end,:)),'b')
%         xlabel('t [s]')
%         ylabel('Alpha [deg]')
%         legend('left','right')
%         title('Angle of attack wingtips')
%         
%         n = n+1;
% 
%         figure(n)
%         hFig = figure(n);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 600 600]);
%         plot(T,radtodeg(alfa_dot_L(end,:)),'r',T,radtodeg(alfa_dot_R(end,:)),'b')
%         xlabel('t [s]')
%         ylabel('Alpha dot [deg/s]')
%         legend('left','right')
%         title('Derivative angle of attack wingtips')
%         
%         n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,FI_vel_b(1,:),'r',T,FA_b(1,:),'b',T,Fg_b(1,:),'g',T,FI_acc_b(1,:),'k',T,ones(N,1)*FA_b_mean(1),'b',T,ones(N,1)*FI_acc_b_mean(1),'k')
        xlabel('t [s]')
        ylabel('F_x')
        subplot(3,1,2); plot(T,FI_vel_b(2,:),'r',T,FA_b(2,:),'b',T,Fg_b(2,:),'g',T,FI_acc_b(2,:),'k',T,ones(N,1)*FA_b_mean(2),'b',T,ones(N,2)*FI_acc_b_mean(2),'k')
        xlabel('t [s]')
        ylabel('F_y')
        subplot(3,1,3); plot(T,FI_vel_b(3,:),'r',T,FA_b(3,:),'b',T,Fg_b(3,:),'g',T,FI_acc_b(3,:),'k',T,ones(N,1)*FA_b_mean(3),'b',T,ones(N,3)*FI_acc_b_mean(3),'k')
        xlabel('t [s]')
        ylabel('F_z')
        legend('FI vel','FA','Fg','FI acc','mean FA','mean FI acc')
        hold off
        [~,h1] = suplabel('Forces (body ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,MI_vel_b(1,:),'r',T,MA_b(1,:),'b',T,MI_acc_b(1,:),'k',T,ones(N,1)*MI_acc_b_mean(1),'k')
        xlabel('t [s]')
        ylabel('M_x')
        subplot(3,1,2); plot(T,MI_vel_b(2,:),'r',T,MA_b(2,:),'b',T,MI_acc_b(2,:),'k',T,ones(N,2)*MI_acc_b_mean(2),'k')
        xlabel('t [s]')
        ylabel('M_y')
        subplot(3,1,3); plot(T,MI_vel_b(3,:),'r',T,MA_b(3,:),'b',T,MI_acc_b(3,:),'k',T,ones(N,3)*MI_acc_b_mean(3),'k')
        xlabel('t [s]')
        ylabel('M_z')
        legend('MI vel','MA','MI acc','mean MI acc')
        hold off
        [~,h1] = suplabel('Moments (body ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,FI_vel_strk(1,:),'r',T,FA_strk(1,:),'b',T,Fg_strk(1,:),'g',T,FI_acc_strk(1,:),'k',T,ones(N,1)*FA_strk_mean(1),'b',T,ones(N,1)*FI_acc_strk_mean(1),'k')
        xlabel('t [s]')
        ylabel('F_x')
        subplot(3,1,2); plot(T,FI_vel_strk(2,:),'r',T,FA_strk(2,:),'b',T,Fg_strk(2,:),'g',T,FI_acc_strk(2,:),'k',T,ones(N,1)*FA_strk_mean(2),'b',T,ones(N,2)*FI_acc_strk_mean(2),'k')
        xlabel('t [s]')
        ylabel('F_y')
        subplot(3,1,3); plot(T,FI_vel_strk(3,:),'r',T,FA_strk(3,:),'b',T,Fg_strk(3,:),'g',T,FI_acc_strk(3,:),'k',T,ones(N,1)*FA_strk_mean(3),'b',T,ones(N,3)*FI_acc_strk_mean(3),'k')
        xlabel('t [s]')
        ylabel('F_z')
        legend('FI vel','FA','Fg','FI acc','mean FA','mean FI acc')
        hold off
        [~,h1] = suplabel('Forces (strokeplane ref. frame)', 't');
        set(h1,'FontSize',10)
        
        n = n+1;

        figure(n)
        hFig = figure(n);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 600 600]);
        hold on
        subplot(3,1,1); plot(T,MI_vel_strk(1,:),'r',T,MA_strk(1,:),'b',T,MI_acc_strk(1,:),'k',T,ones(N,1)*MI_acc_strk_mean(1),'k')
        xlabel('t [s]')
        ylabel('M_x')
        subplot(3,1,2); plot(T,MI_vel_strk(2,:),'r',T,MA_strk(2,:),'b',T,MI_acc_strk(2,:),'k',T,ones(N,2)*MI_acc_strk_mean(2),'k')
        xlabel('t [s]')
        ylabel('M_y')
        subplot(3,1,3); plot(T,MI_vel_strk(3,:),'r',T,MA_strk(3,:),'b',T,MI_acc_strk(3,:),'k',T,ones(N,3)*MI_acc_strk_mean(3),'k')
        xlabel('t [s]')
        ylabel('M_z')
        legend('MI vel','MA','MI acc','mean MI acc')
        hold off
        [~,h1] = suplabel('Moments (strokeplane ref. frame)', 't');
        set(h1,'FontSize',10)
        
    end

end

