function Sphere_plot_raw_filt2(settings,seq_nr,qL1_r,qL2_r,qL3_r,qL4_r,qR1_r,qR2_r,qR3_r,qR4_r,qL1_f2,qL2_f2,qL3_f2,qL4_f2,qR1_f2,qR2_f2,qR3_f2,qR4_f2,wing_l,fig_nr1,save_on_off)
    
    Lwingtip_r = zeros(length(qL1_r),3);

    Lwingtip_f2 = zeros(length(qL1_r),3);
    
    Lwt = wing_l.*[0; -1; 0];
    
    Rwingtip_r = zeros(length(qR1_r),3);

    Rwingtip_f2 = zeros(length(qR1_r),3);
    
    Rwt = wing_l.*[0; 1; 0];
        
    % Calculate wing_tip path
    
    for j = 1:length(qL1_r)
        
        DCM_L_r = quat2matNEW([qL1_r(j) qL2_r(j) qL3_r(j) qL4_r(j)]);

        Lwingtip_r(j,:) = DCM_L_r*Lwt;
        
        DCM_R_r = quat2matNEW([qR1_r(j) qR2_r(j) qR3_r(j) qR4_r(j)]);

        Rwingtip_r(j,:) = DCM_R_r*Rwt;
        
    end

    for j = 1:length(qL1_r)
        
        DCM_L_f2 = quat2matNEW([qL1_f2(j) qL2_f2(j) qL3_f2(j) qL4_f2(j)]);

        Lwingtip_f2(j,:) = DCM_L_f2*Lwt;
        
        DCM_R_f2 = quat2matNEW([qR1_f2(j) qR2_f2(j) qR3_f2(j) qR4_f2(j)]);

        Rwingtip_f2(j,:) = DCM_R_f2*Rwt;
        
    end
    
       
    k = 5;
    n = 2^k-1;
    [x,y,z] = sphere(n);
    
%     figure(fig_nr1)
%     surf(wing_l.*x,wing_l.*y,wing_l.*z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold on
%     
%     figure(fig_nr1)
%     plot3(Lwingtip_r(:,1),Lwingtip_r(:,2),Lwingtip_r(:,3),'r')
%     axis equal
%     hold on
%     
%     figure(fig_nr1)
%     plot3(Rwingtip_r(:,1),Rwingtip_r(:,2),Rwingtip_r(:,3),'b')
%     axis equal
%     title('Sphere plot wingtip-path unfiltered left and right wing')
%     xlabel('x body [mm]')
%     ylabel('y body [mm]')
%     zlabel('z body [mm]')
%     hold off
%     
%     
%     
%     figure(fig_nr2)
%     surf(wing_l.*x,wing_l.*y,wing_l.*z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold on
%     
%     figure(fig_nr2)
%     plot3(Lwingtip_f2(:,1),Lwingtip_f2(:,2),Lwingtip_f2(:,3),'r')
%     axis equal
%     hold on
%     
%     figure(fig_nr2)
%     plot3(Rwingtip_f2(:,1),Rwingtip_f2(:,2),Rwingtip_f2(:,3),'b')
%     axis equal
%     title('Sphere plot wingtip-path filtered left and right wing')
%     xlabel('x body [mm]')
%     ylabel('y body [mm]')
%     zlabel('z body [mm]')
%     hold off
    
    
    
    
    figure(fig_nr1)
    surf(wing_l.*x,wing_l.*y,wing_l.*z,'FaceColor','black','EdgeColor','none');
    alpha(0.2)
    axis equal
    hold on
    
    figure(fig_nr1)
    plot3(Lwingtip_r(:,1),Lwingtip_r(:,2),Lwingtip_r(:,3),'g')
    axis equal
    hold on
    
    figure(fig_nr1)
    plot3(Rwingtip_r(:,1),Rwingtip_r(:,2),Rwingtip_r(:,3),'g')
    axis equal
    hold on
    
    figure(fig_nr1)
    plot3(Lwingtip_f2(:,1),Lwingtip_f2(:,2),Lwingtip_f2(:,3),'r')
    axis equal
    hold on
    
    figure(fig_nr1)
    plot3(Rwingtip_f2(:,1),Rwingtip_f2(:,2),Rwingtip_f2(:,3),'r')
    axis equal
    title('Sphere plot wingtip-path of filtered and unfiltered left and right wing')
    xlabel('x body [mm]')
    ylabel('y body [mm]')
    zlabel('z body [mm]')
    hold off
    
    % Save plots
    
    if save_on_off == 1
    
    saveas(fig_nr1, [char(settings.plot_folders(1)) '/' char(settings.sequence_names(seq_nr)) '/sphere1'], 'fig')
    
    saveas(fig_nr1, [char(settings.plot_folders(2)) '/sphere_plot/sphere1_' int2str(seq_nr)], 'fig')
    
    end


end

