function Kalman_filt_plots(settings,pathDB)

    seq_nr = 11;
    
    start = find(isnan(pathDB.filt.xyz(:,1,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.filt.xyz(:,1,seq_nr))==0, 1, 'last' );
    
    start2 = start;
    stop2 = start2+300;
    
    t2 = pathDB.t(start2:stop2);
    t = pathDB.t(start:stop);
    dt = pathDB.dt;
    
    xyz_raw = pathDB.raw.xyz(start2:stop2,:,seq_nr);
    xyz_filt = pathDB.filt.xyz(start2:stop2,:,seq_nr);
    
    x_min = min([xyz_raw(:,1); xyz_filt(:,1)]);
    x_max = max([xyz_raw(:,1); xyz_filt(:,1)]);
    y_min = min([xyz_raw(:,2); xyz_filt(:,2)]);
    y_max = max([xyz_raw(:,2); xyz_filt(:,2)]);
    z_min = min([xyz_raw(:,3); xyz_filt(:,3)]);
    z_max = max([xyz_raw(:,3); xyz_filt(:,3)]);
    
    v_filt = pathDB.filt.uvw(start:stop,:,seq_nr);
    
    a_filt = pathDB.filt.a_xyz(start:stop,:,seq_nr);
    
    qb_raw   = pathDB.raw.qB(start:stop,:,seq_nr);
    
    qb_filt  = pathDB.filt.qB(start:stop,:,seq_nr);
    
    wb_filt = pathDB.filt.wB(start:stop,:,seq_nr);
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1400 800]);
    hold on
    subplot(3,2,[1; 3; 5]);     plot3(xyz_raw(:,1),xyz_raw(:,2),xyz_raw(:,3))
    title('Trajectory')
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    xlim([x_min x_max])
    ylim([y_min y_max])
    zlim([z_min z_max])
    grid on
    axis equal
    subplot(3,2,2);             plot(t2,xyz_raw(:,1))
    title('x coordinates')
    xlabel('t [s]')
    ylabel('x [mm]')
    subplot(3,2,4);             plot(t2,xyz_raw(:,2))
    title('y coordinates')
    xlabel('t [s]')
    ylabel('y [mm]')
    subplot(3,2,6);             plot(t2,xyz_raw(:,3))
    title('z coordinates')
    xlabel('t [s]')
    ylabel('z [mm]')
    hold off
    [~,h1] = suplabel('Unfiltered position body c.g.', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/xyz_unfilt'],'fig')
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FX_forward.eps'] ,'-depsc2');

    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1400 800]);
    hold on
    subplot(3,2,[1; 3; 5]); hold on
    plot3(xyz_filt(:,1),xyz_filt(:,2),xyz_filt(:,3),'r')
    plot3(xyz_raw(:,1),xyz_raw(:,2),xyz_raw(:,3),'b')
    hold off
    title('Trajectory')
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    xlim([x_min x_max])
    ylim([y_min y_max])
    zlim([z_min z_max])
    legend('filtered','unfiltered')
    grid on
    axis equal
    subplot(3,2,2); hold on
    plot(t2,xyz_raw(:,1),'b')
    plot(t2,xyz_filt(:,1),'r')
    hold off
    title('x coordinates')
    xlabel('t [s]')
    ylabel('x [mm]')
    subplot(3,2,4); hold on
    plot(t2,xyz_raw(:,2),'b')
    plot(t2,xyz_filt(:,2),'r')
    hold off
    title('y coordinates')
    xlabel('t [s]')
    ylabel('y [mm]')
    subplot(3,2,6); hold on
    plot(t2,xyz_raw(:,3),'b')
    plot(t2,xyz_filt(:,3),'r')
    hold off
    title('z coordinates')
    xlabel('t [s]')
    ylabel('z [mm]')
    hold off
    [~,h1] = suplabel('Unfiltered and filtered position body c.g.', 't');
    set(h1,'FontSize',10)
%     print ([char(settings.plot_loc) '/Maneuvering_wing_kin/FX_forward.eps'] ,'-depsc2');
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/xyz_filt'],'fig')

    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 600]);
    hold on
    subplot(3,1,1); plot(t,v_filt(:,1),'r')
    ylabel('v_x [mm/s]')
    subplot(3,1,2); plot(t,v_filt(:,2),'r')
    ylabel('v_y [mm/s]')
    subplot(3,1,3); plot(t,v_filt(:,3),'r')
    xlabel('t [s]')
    ylabel('v_z [mm/s]')
    hold off
    [~,h1] = suplabel('Filtered velocity body c.g.', 't');
    set(h1,'FontSize',10)    
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/v_filt'],'fig')
    
    figure(4)
    hFig = figure(4);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 600]);
    hold on
    subplot(3,1,1); plot(t,a_filt(:,1),'r')
    ylabel('a_x [mm/s^2]')
    subplot(3,1,2); plot(t,a_filt(:,2),'r')
    ylabel('a_y [mm/s^2]')
    subplot(3,1,3); plot(t,a_filt(:,3),'r')
    xlabel('t [s]')
    ylabel('a_z [mm/s^2]')
    hold off
    [~,h1] = suplabel('Filtered acceleration body c.g.', 't');
    set(h1,'FontSize',10)    
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/a_filt'],'fig')
    
    figure(5)
    hFig = figure(5);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 600]);
    hold on
    subplot(4,1,1); plot(t,qb_raw(:,1),'b',t,qb_filt(:,1),'r')
    ylabel('q_{b1}')
    subplot(4,1,2); plot(t,qb_raw(:,2),'b',t,qb_filt(:,2),'r')
    ylabel('q_{b2}')
    subplot(4,1,3); plot(t,qb_raw(:,3),'b',t,qb_filt(:,3),'r')
    ylabel('q_{b3}')
    subplot(4,1,4); plot(t,qb_raw(:,4),'b',t,qb_filt(:,4),'r')
    xlabel('t [s]')
    ylabel('q_{b4}')
    hold off
    [~,h1] = suplabel('Unfiltered and filtered body quaternion', 't');
    set(h1,'FontSize',10)    
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/qb_filt'],'fig')

    figure(6)
    hFig = figure(6);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 600]);
    hold on
    subplot(3,1,1); plot(t,wb_filt(:,1),'r')
    ylabel('\omega_{b1} [rad/s]')
    subplot(3,1,2); plot(t,wb_filt(:,2),'r')
    ylabel('\omega_{b2} [rad/s]')
    subplot(3,1,3); plot(t,wb_filt(:,3),'r')
    xlabel('t [s]')
    ylabel('\omega_{b3} [rad/s]')
    hold off
    [~,h1] = suplabel('Filtered body angular velocity', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/wb_filt'],'fig')
    
    R_strk = pathDB.rot_mat.Rstr;
    
    Rb = R_strk';
    
    N = stop-start+1;
    
    RL_raw = zeros(3,3,N);
    RR_raw = zeros(3,3,N);
    
    RL_filt = zeros(3,3,N);
    RR_filt = zeros(3,3,N);
    
    qL_raw = pathDB.raw.qL(start:stop,:,seq_nr);
    qR_raw = pathDB.raw.qR(start:stop,:,seq_nr);
    
    qL_filt = pathDB.filt.qL(start:stop,:,seq_nr);
    qR_filt = pathDB.filt.qR(start:stop,:,seq_nr);
    
    for i = 1:N
        
        RL_raw(:,:,i)   = quat2mat(qL_raw(i,:));
        RR_raw(:,:,i)   = quat2mat(qR_raw(i,:));
        RL_filt(:,:,i)  = quat2mat(qL_filt(i,:));
        RR_filt(:,:,i)  = quat2mat(qR_filt(i,:));
        
    end
    
    figure(7)
    hFig = figure(7);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 1000]);
    hold on
    sphere_plot_3D( settings, pathDB, Rb, RL_raw, RR_raw, RL_filt, RR_filt, seq_nr)
    title('Sphere plot wingtip trajectory')
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    axis equal
    grid on
    hold off
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/sphere_plot'],'fig')
    
    figure(8)
    hFig = figure(8);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 600]);
    hold on
    subplot(4,1,1); plot(t2,qL_raw(1:(stop2-start2+1),1),'b',t2,qL_filt(1:(stop2-start2+1),1),'r')
    ylabel('q_{L1}')
    subplot(4,1,2); plot(t2,qL_raw(1:(stop2-start2+1),2),'b',t2,qL_filt(1:(stop2-start2+1),2),'r')
    ylabel('q_{L2}')
    subplot(4,1,3); plot(t2,qL_raw(1:(stop2-start2+1),3),'b',t2,qL_filt(1:(stop2-start2+1),3),'r')
    ylabel('q_{L3}')
    subplot(4,1,4); plot(t2,qL_raw(1:(stop2-start2+1),4),'b',t2,qL_filt(1:(stop2-start2+1),4),'r')
    xlabel('t [s]')
    ylabel('q_{L4}')
    hold off
    [~,h1] = suplabel('Unfiltered and filtered left wing quaternion', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/qL_filt'],'fig')
    
    wL_filt = pathDB.filt.wL(start:stop,:,seq_nr);
    
    figure(9)
    hFig = figure(9);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 600]);
    hold on
    subplot(3,1,1); plot(t2,wL_filt(1:(stop2-start2+1),1),'r')
    ylabel('\omega_{L1} [rad/s]')
    subplot(3,1,2); plot(t2,wL_filt(1:(stop2-start2+1),2),'r')
    ylabel('\omega_{L2} [rad/s]')
    subplot(3,1,3); plot(t2,wL_filt(1:(stop2-start2+1),3),'r')
    xlabel('t [s]')
    ylabel('\omega_{L3} [rad/s]')
    hold off
    [~,h1] = suplabel('Filtered angular velocity left wing', 't');
    set(h1,'FontSize',10)
    saveas(hFig,[char(settings.plot_loc) '/Kalman_filt_plots/wL_filt'],'fig')
    
end