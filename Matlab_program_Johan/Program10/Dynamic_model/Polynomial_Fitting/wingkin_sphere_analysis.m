function wingkin_sphere_analysis(settings,pathDB)


    dt = pathDB.t(2)-pathDB.t(1);
    
    nr_of_seq = size(pathDB.x,2);
    
    LWT = {};
    
    RWT = {};
    
    LWT.A = 0;
    LWT.B = 0;
    LWT.C = 0;
    LWT.D = 0;
    
    RWT.A = 0;
    RWT.B = 0;
    RWT.C = 0;
    RWT.D = 0;

    
    STRK = {};
    
    STRK.w_x = 0;
    STRK.w_y = 0;
    STRK.w_z = 0;
    
    STRK.w_dot_x = 0;
    STRK.w_dot_y = 0;
    STRK.w_dot_z = 0;
    
    STRK.Vn = 0;
    STRK.Vt = 0;
    
    STRK.An = 0;
    STRK.At = 0;
        


    for i =1:nr_of_seq
    
    seq_nr = i;
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );

    % Analyze the wing kinematics on the spherical plane:
    
    down_time = pathDB.down_time_L(:,:,seq_nr);
    up_time = pathDB.up_time_L(:,:,seq_nr);
    
    nr_wb = find(isnan(down_time(:,1))==0, 1, 'last' );
    
    theta_L = pathDB.theta_L(start:stop,seq_nr);
    eta_L = pathDB.eta_L(start:stop,seq_nr);
    phi_L = pathDB.phi_L(start:stop,seq_nr);
    
    theta_R = pathDB.theta_R(start:stop,seq_nr);
    eta_R = pathDB.eta_R(start:stop,seq_nr);
    phi_R = pathDB.phi_R(start:stop,seq_nr);
% 
%     qb1_ = pathDB.qb1_filt(start:stop,seq_nr);
%     qb2 = pathDB.qb2_filt(start:stop,seq_nr);
%     qb3 = pathDB.qb3_filt(start:stop,seq_nr);
%     qb4 = pathDB.qb4_filt(start:stop,seq_nr);
    
    qL1 = pathDB.qL1_filt2(start:stop,seq_nr);
    qL2 = pathDB.qL2_filt2(start:stop,seq_nr);
    qL3 = pathDB.qL3_filt2(start:stop,seq_nr);
    qL4 = pathDB.qL4_filt2(start:stop,seq_nr);
    
    qR1 = pathDB.qR1_filt2(start:stop,seq_nr);
    qR2 = pathDB.qR2_filt2(start:stop,seq_nr);
    qR3 = pathDB.qR3_filt2(start:stop,seq_nr);
    qR4 = pathDB.qR4_filt2(start:stop,seq_nr);
    
    
    A_loc_L =  zeros(4,nr_wb);
    B_loc_L =  zeros(4,nr_wb);
    C_loc_L =  zeros(4,nr_wb);
    D_loc_L =  zeros(4,nr_wb);
    
    A_loc_R =  zeros(4,nr_wb);
    B_loc_R =  zeros(4,nr_wb);
    C_loc_R =  zeros(4,nr_wb);
    D_loc_R =  zeros(4,nr_wb);    
    
    for j = 1:nr_wb
        
        down_time_end = find(isnan(down_time(j,:))==0, 1, 'last' );
       
        up_time_end =  find(isnan(up_time(j,:))==0, 1, 'last' );
        
        theta_L_down = theta_L(down_time(j,1:down_time_end));
       
        theta_R_down = theta_R(down_time(j,1:down_time_end));
       
        theta_L_up = theta_L(up_time(j,1:up_time_end));
       
        theta_R_up = theta_R(up_time(j,1:up_time_end));

        
        A_loc_id_L = down_time(j,1);
        
        [~, B_loc_id_L] = min(theta_L_down);
        
        C_loc_id_L = up_time(j,1);
        
        [~, D_loc_id_L] = min(theta_L_up);
        
        A_loc_id_R = down_time(j,1);
        
        [~, B_loc_id_R] = min(theta_R_down);
        
        C_loc_id_R = up_time(j,1);
        
        [~, D_loc_id_R] = min(theta_R_up);

%         [~, A_loc_id_L] = max(theta_L_down);
%         
%         [~, B_loc_id_L] = min(theta_L_down);
%         
%         [~, C_loc_id_L] = max(theta_L_up);
%         
%         [~, D_loc_id_L] = min(theta_L_up);
%         
%         [~, A_loc_id_R] = max(theta_R_down);
%         
%         [~, B_loc_id_R] = min(theta_R_down);
%         
%         [~, C_loc_id_R] = max(theta_R_up);
%         
%         [~, D_loc_id_R] = min(theta_R_up);
        
        A_loc_L(:,j) = [qL1(A_loc_id_L); qL2(A_loc_id_L); qL3(A_loc_id_L); qL4(A_loc_id_L)];
        B_loc_L(:,j) = [qL1(down_time(j,1)-1+B_loc_id_L); qL2(down_time(j,1)-1+B_loc_id_L); qL3(down_time(j,1)-1+B_loc_id_L); qL4(down_time(j,1)-1+B_loc_id_L)];
        C_loc_L(:,j) = [qL1(C_loc_id_L); qL2(C_loc_id_L); qL3(C_loc_id_L); qL4(C_loc_id_L)];
        D_loc_L(:,j) = [qL1(up_time(j,1)-1+D_loc_id_L); qL2(up_time(j,1)-1+D_loc_id_L); qL3(up_time(j,1)-1+D_loc_id_L); qL4(up_time(j,1)-1+D_loc_id_L)];

        A_loc_R(:,j) = [qR1(A_loc_id_R); qR2(A_loc_id_R); qR3(A_loc_id_R); qR4(A_loc_id_R)];
        B_loc_R(:,j) = [qR1(down_time(j,1)-1+B_loc_id_R); qR2(down_time(j,1)-1+B_loc_id_R); qR3(down_time(j,1)-1+B_loc_id_R); qR4(down_time(j,1)-1+B_loc_id_R)];
        C_loc_R(:,j) = [qR1(C_loc_id_R); qR2(C_loc_id_R); qR3(C_loc_id_R); qR4(C_loc_id_R)];
        D_loc_R(:,j) = [qR1(up_time(j,1)-1+D_loc_id_R); qR2(up_time(j,1)-1+D_loc_id_R); qR3(up_time(j,1)-1+D_loc_id_R); qR4(up_time(j,1)-1+D_loc_id_R)];


%         A_loc_L(:,j) = [qL1(down_time(j,1)-1+A_loc_id_L); qL2(down_time(j,1)-1+A_loc_id_L); qL3(down_time(j,1)-1+A_loc_id_L); qL4(down_time(j,1)-1+A_loc_id_L)];
%         B_loc_L(:,j) = [qL1(down_time(j,1)-1+B_loc_id_L); qL2(down_time(j,1)-1+B_loc_id_L); qL3(down_time(j,1)-1+B_loc_id_L); qL4(down_time(j,1)-1+B_loc_id_L)];
%         C_loc_L(:,j) = [qL1(up_time(j,1)-1+C_loc_id_L); qL2(up_time(j,1)-1+C_loc_id_L); qL3(up_time(j,1)-1+C_loc_id_L); qL4(up_time(j,1)-1+C_loc_id_L)];
%         D_loc_L(:,j) = [qL1(up_time(j,1)-1+D_loc_id_L); qL2(up_time(j,1)-1+D_loc_id_L); qL3(up_time(j,1)-1+D_loc_id_L); qL4(up_time(j,1)-1+D_loc_id_L)];
% 
%         A_loc_R(:,j) = [qR1(down_time(j,1)-1+A_loc_id_R); qR2(down_time(j,1)-1+A_loc_id_R); qR3(down_time(j,1)-1+A_loc_id_R); qR4(down_time(j,1)-1+A_loc_id_R)];
%         B_loc_R(:,j) = [qR1(down_time(j,1)-1+B_loc_id_R); qR2(down_time(j,1)-1+B_loc_id_R); qR3(down_time(j,1)-1+B_loc_id_R); qR4(down_time(j,1)-1+B_loc_id_R)];
%         C_loc_R(:,j) = [qR1(up_time(j,1)-1+C_loc_id_R); qR2(up_time(j,1)-1+C_loc_id_R); qR3(up_time(j,1)-1+C_loc_id_R); qR4(up_time(j,1)-1+C_loc_id_R)];
%         D_loc_R(:,j) = [qR1(up_time(j,1)-1+D_loc_id_R); qR2(up_time(j,1)-1+D_loc_id_R); qR3(up_time(j,1)-1+D_loc_id_R); qR4(up_time(j,1)-1+D_loc_id_R)];
        
    end
    

        
    Lwt = [0; -1; 0];
    Rwt = [0; 1; 0];
    
    Lwingtip_A = zeros(nr_wb,3);
    Rwingtip_A = zeros(nr_wb,3);
    Lwingtip_B = zeros(nr_wb,3);
    Rwingtip_B = zeros(nr_wb,3);
    Lwingtip_C = zeros(nr_wb,3);
    Rwingtip_C = zeros(nr_wb,3);
    Lwingtip_D = zeros(nr_wb,3);
    Rwingtip_D = zeros(nr_wb,3);

    for j = 1:nr_wb
        
        DCM_L_A = quat2matNEW(A_loc_L(:,j));
        DCM_R_A = quat2matNEW(A_loc_R(:,j));
        DCM_L_B = quat2matNEW(B_loc_L(:,j));
        DCM_R_B = quat2matNEW(B_loc_R(:,j));
        DCM_L_C = quat2matNEW(C_loc_L(:,j));
        DCM_R_C = quat2matNEW(C_loc_R(:,j));
        DCM_L_D = quat2matNEW(D_loc_L(:,j));
        DCM_R_D = quat2matNEW(D_loc_R(:,j));
        

        Lwingtip_A(j,:) = DCM_L_A*Lwt;       
        Rwingtip_A(j,:) = DCM_R_A*Rwt;
        Lwingtip_B(j,:) = DCM_L_B*Lwt;       
        Rwingtip_B(j,:) = DCM_R_B*Rwt;
        Lwingtip_C(j,:) = DCM_L_C*Lwt;       
        Rwingtip_C(j,:) = DCM_R_C*Rwt;
        Lwingtip_D(j,:) = DCM_L_D*Lwt;       
        Rwingtip_D(j,:) = DCM_R_D*Rwt;
       
    end
    
    
    % Strokeplane dynamics:

    [ stroke_var ] = Strokeplane_dynamics( settings, pathDB, i );
        
    if i == 1
       
       LWT.A = Lwingtip_A;
       LWT.B = Lwingtip_B;
       LWT.C = Lwingtip_C;
       LWT.D = Lwingtip_D;
       
       RWT.A = Rwingtip_A;
       RWT.B = Rwingtip_B;
       RWT.C = Rwingtip_C;
       RWT.D = Rwingtip_D;
       
       STRK.w_x = stroke_var.Omega_strk(1,:)';
       STRK.w_y = stroke_var.Omega_strk(2,:)';
       STRK.w_z = stroke_var.Omega_strk(3,:)';

       STRK.w_dot_x = stroke_var.Omega_dot_strk(:,1)';
       STRK.w_dot_y = stroke_var.Omega_dot_strk(:,2)';
       STRK.w_dot_z = stroke_var.Omega_dot_strk(:,3)';

       STRK.Vn = stroke_var.Vn';
       STRK.Vt = stroke_var.Vt';

       STRK.An = stroke_var.An';
       STRK.At = stroke_var.At';
        
    else
       
       LWT.A = [LWT.A; Lwingtip_A];
       LWT.B = [LWT.B; Lwingtip_B];
       LWT.C = [LWT.C; Lwingtip_C];
       LWT.D = [LWT.D; Lwingtip_D];
        
       RWT.A = [RWT.A; Rwingtip_A];
       RWT.B = [RWT.B; Rwingtip_B];
       RWT.C = [RWT.C; Rwingtip_C];
       RWT.D = [RWT.D; Rwingtip_D];
       
       STRK.w_x = [STRK.w_x; stroke_var.Omega_strk(1,:)'];
       STRK.w_y = [STRK.w_y; stroke_var.Omega_strk(2,:)'];
       STRK.w_z = [STRK.w_z; stroke_var.Omega_strk(3,:)'];

       STRK.w_dot_x = [STRK.w_dot_x; stroke_var.Omega_dot_strk(:,1)'];
       STRK.w_dot_y = [STRK.w_dot_y; stroke_var.Omega_dot_strk(:,2)'];
       STRK.w_dot_z = [STRK.w_dot_z; stroke_var.Omega_dot_strk(:,3)'];

       STRK.Vn = [STRK.Vn; stroke_var.Vn'];
       STRK.Vt = [STRK.Vt; stroke_var.Vt'];

       STRK.An = [STRK.An; stroke_var.An'];
       STRK.At = [STRK.At; stroke_var.At'];
               
    end

       
  
    
    
    k = 5;
    n = 2^k-1;
    [x,y,z] = sphere(n);   
    
    
    C_map = jet(200);
    
   
    
    figure()
    hold on
        for k = 1:nr_wb
            if stroke_var.Omega_dot_strk(1,k) == 0
                color_code_wx = 100;
            elseif stroke_var.Omega_dot_strk(1,k) > 0
                color_code_wx = 101+round(stroke_var.Omega_dot_strk(1,k)/max(stroke_var.Omega_dot_strk(1,:))*99);
            elseif  stroke_var.Omega_dot_strk(1,k) < 0
                color_code_wx = 1+round(stroke_var.Omega_dot_strk(1,k)/min(stroke_var.Omega_dot_strk(1,:))*99);
            end
            plot3(LWT.C(k,1),LWT.C(k,2),LWT.C(k,3),'o','Color',C_map(color_code_wx,:))
            plot3(RWT.C(k,1),RWT.C(k,2),RWT.C(k,3),'o','Color',C_map(color_code_wx,:))

        end
    surf(x,y,z,'FaceColor','black','EdgeColor','none');
    alpha(0.2)
    axis equal
    colorbar
    hold off
    
    
    
    
%     figure(seq_nr)
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Lwingtip_A(k,1),Lwingtip_A(k,2),Lwingtip_A(k,3),'o','Color','r')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Rwingtip_A(k,1),Rwingtip_A(k,2),Rwingtip_A(k,3),'o','Color','r')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Lwingtip_B(k,1),Lwingtip_B(k,2),Lwingtip_B(k,3),'o','Color','b')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Rwingtip_B(k,1),Rwingtip_B(k,2),Rwingtip_B(k,3),'o','Color','b')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Lwingtip_C(k,1),Lwingtip_C(k,2),Lwingtip_C(k,3),'o','Color','g')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Rwingtip_C(k,1),Rwingtip_C(k,2),Rwingtip_C(k,3),'o','Color','g')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Lwingtip_D(k,1),Lwingtip_D(k,2),Lwingtip_D(k,3),'o','Color','m')
%     end
%     axis equal
%     hold on
%     
%     figure(seq_nr)
%     for k = 1:nr_wb
%     plot3(Rwingtip_D(k,1),Rwingtip_D(k,2),Rwingtip_D(k,3),'o','Color','m')
%     end
%     axis equal
%     title('Sphere plot wingtip-path of filtered and unfiltered left and right wing')
%     xlabel('x body [mm]')
%     ylabel('y body [mm]')
%     zlabel('z body [mm]')
%     hold off

    pause
    
    end
 
%     
%     color_code_w_x = round(w_x_sort./max(w_x_sort));
    
%     plot(xi,stroke_wb_L_interp,'-','color',cmap_jet(color_code_w_x,:),'linewidth',linewidth_timelines)

%     [w_x_sort, w_x_sort_id] = sort(STRK.w_x);
%     [w_y_sort, w_y_sort_id] = sort(STRK.w_y);
%     [w_z_sort, w_z_sort_id] = sort(STRK.w_z);
%     
%     [w_dot_x_sort, w_dot_x_sort_id] = sort(STRK.w_dot_x);
%     [w_dot_y_sort, w_dot_y_sort_id] = sort(STRK.w_dot_y);
%     [w_dot_z_sort, w_dot_z_sort_id] = sort(STRK.w_dot_z);
%     
%     [Vn_sort, Vn_sort_id] = sort(STRK.Vn);
%     [Vt_sort, Vt_sort_id] = sort(STRK.Vt);
%     
%     [An_sort, An_sort_id] = sort(STRK.An);
%     [At_sort, At_sort_id] = sort(STRK.At);
% 
% 
%     
%     
%     k = 5;
%     n = 2^k-1;
%     [x,y,z] = sphere(n);   
%     
%     
%     C_map = jet(200);
%     
%    
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_x((j*(k-1))+k,1) == 0
%                 color_code_wx = 100;
%             elseif STRK.w_x((j*(k-1))+k,1) > 0
%                 color_code_wx = 101+round((STRK.w_x((j*(k-1))+k,1)/max(STRK.w_x))*99);
%             elseif STRK.w_x((j*(k-1))+k,1) < 0
%                 color_code_wx = 1+round((STRK.w_x((j*(k-1))+k,1)/min(STRK.w_x))*99);
%             end
%             plot3(LWT.A((j*(k-1))+k,1),LWT.A((j*(k-1))+k,2),LWT.A((j*(k-1))+k,3),'o','Color',C_map(color_code_wx,:))
%             plot3(RWT.A((j*(k-1))+k,1),RWT.A((j*(k-1))+k,2),RWT.A((j*(k-1))+k,3),'o','Color',C_map(color_code_wx,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_y((j*(k-1))+k,1) == 0
%                 color_code_wy = 100;
%             elseif STRK.w_y((j*(k-1))+k,1) > 0
%                 color_code_wy = 101+round((STRK.w_y((j*(k-1))+k,1)/max(STRK.w_y))*99);
%             elseif STRK.w_y((j*(k-1))+k,1) < 0
%                 color_code_wy = 1+round((STRK.w_y((j*(k-1))+k,1)/min(STRK.w_y))*99);
%             end
%             plot3(LWT.A((j*(k-1))+k,1),LWT.A((j*(k-1))+k,2),LWT.A((j*(k-1))+k,3),'o','Color',C_map(color_code_wy,:))
%             plot3(RWT.A((j*(k-1))+k,1),RWT.A((j*(k-1))+k,2),RWT.A((j*(k-1))+k,3),'o','Color',C_map(color_code_wy,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_z((j*(k-1))+k,1) == 0
%                 color_code_wz = 100;
%             elseif STRK.w_z((j*(k-1))+k,1) > 0
%                 color_code_wz = 101+round((STRK.w_z((j*(k-1))+k,1)/max(STRK.w_z))*99);
%             elseif STRK.w_z((j*(k-1))+k,1) < 0
%                 color_code_wz = 1+round((STRK.w_z((j*(k-1))+k,1)/min(STRK.w_z))*99);
%             end
%             plot3(LWT.A((j*(k-1))+k,1),LWT.A((j*(k-1))+k,2),LWT.A((j*(k-1))+k,3),'o','Color',C_map(color_code_wz,:))
%             plot3(RWT.A((j*(k-1))+k,1),RWT.A((j*(k-1))+k,2),RWT.A((j*(k-1))+k,3),'o','Color',C_map(color_code_wz,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     
%         figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_x((j*(k-1))+k,1) == 0
%                 color_code_wx = 100;
%             elseif STRK.w_x((j*(k-1))+k,1) > 0
%                 color_code_wx = 101+round((STRK.w_x((j*(k-1))+k,1)/max(STRK.w_x))*99);
%             elseif STRK.w_x((j*(k-1))+k,1) < 0
%                 color_code_wx = 1+round((STRK.w_x((j*(k-1))+k,1)/min(STRK.w_x))*99);
%             end
%             plot3(LWT.B((j*(k-1))+k,1),LWT.B((j*(k-1))+k,2),LWT.B((j*(k-1))+k,3),'o','Color',C_map(color_code_wx,:))
%             plot3(RWT.B((j*(k-1))+k,1),RWT.B((j*(k-1))+k,2),RWT.B((j*(k-1))+k,3),'o','Color',C_map(color_code_wx,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_y((j*(k-1))+k,1) == 0
%                 color_code_wy = 100;
%             elseif STRK.w_y((j*(k-1))+k,1) > 0
%                 color_code_wy = 101+round((STRK.w_y((j*(k-1))+k,1)/max(STRK.w_y))*99);
%             elseif STRK.w_y((j*(k-1))+k,1) < 0
%                 color_code_wy = 1+round((STRK.w_y((j*(k-1))+k,1)/min(STRK.w_y))*99);
%             end
%             plot3(LWT.B((j*(k-1))+k,1),LWT.B((j*(k-1))+k,2),LWT.B((j*(k-1))+k,3),'o','Color',C_map(color_code_wy,:))
%             plot3(RWT.B((j*(k-1))+k,1),RWT.B((j*(k-1))+k,2),RWT.B((j*(k-1))+k,3),'o','Color',C_map(color_code_wy,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_z((j*(k-1))+k,1) == 0
%                 color_code_wz = 100;
%             elseif STRK.w_z((j*(k-1))+k,1) > 0
%                 color_code_wz = 101+round((STRK.w_z((j*(k-1))+k,1)/max(STRK.w_z))*99);
%             elseif STRK.w_z((j*(k-1))+k,1) < 0
%                 color_code_wz = 1+round((STRK.w_z((j*(k-1))+k,1)/min(STRK.w_z))*99);
%             end
%             plot3(LWT.B((j*(k-1))+k,1),LWT.B((j*(k-1))+k,2),LWT.B((j*(k-1))+k,3),'o','Color',C_map(color_code_wz,:))
%             plot3(RWT.B((j*(k-1))+k,1),RWT.B((j*(k-1))+k,2),RWT.B((j*(k-1))+k,3),'o','Color',C_map(color_code_wz,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_x((j*(k-1))+k,1) == 0
%                 color_code_wx = 100;
%             elseif STRK.w_x((j*(k-1))+k,1) > 0
%                 color_code_wx = 101+round((STRK.w_x((j*(k-1))+k,1)/max(STRK.w_x))*99);
%             elseif STRK.w_x((j*(k-1))+k,1) < 0
%                 color_code_wx = 1+round((STRK.w_x((j*(k-1))+k,1)/min(STRK.w_x))*99);
%             end
%             plot3(LWT.C((j*(k-1))+k,1),LWT.C((j*(k-1))+k,2),LWT.C((j*(k-1))+k,3),'o','Color',C_map(color_code_wx,:))
%             plot3(RWT.C((j*(k-1))+k,1),RWT.C((j*(k-1))+k,2),RWT.C((j*(k-1))+k,3),'o','Color',C_map(color_code_wx,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_y((j*(k-1))+k,1) == 0
%                 color_code_wy = 100;
%             elseif STRK.w_y((j*(k-1))+k,1) > 0
%                 color_code_wy = 101+round((STRK.w_y((j*(k-1))+k,1)/max(STRK.w_y))*99);
%             elseif STRK.w_y((j*(k-1))+k,1) < 0
%                 color_code_wy = 1+round((STRK.w_y((j*(k-1))+k,1)/min(STRK.w_y))*99);
%             end
%             plot3(LWT.C((j*(k-1))+k,1),LWT.C((j*(k-1))+k,2),LWT.C((j*(k-1))+k,3),'o','Color',C_map(color_code_wy,:))
%             plot3(RWT.C((j*(k-1))+k,1),RWT.C((j*(k-1))+k,2),RWT.C((j*(k-1))+k,3),'o','Color',C_map(color_code_wy,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             if STRK.w_z((j*(k-1))+k,1) == 0
%                 color_code_wz = 100;
%             elseif STRK.w_z((j*(k-1))+k,1) > 0
%                 color_code_wz = 101+round((STRK.w_z((j*(k-1))+k,1)/max(STRK.w_z))*99);
%             elseif STRK.w_z((j*(k-1))+k,1) < 0
%                 color_code_wz = 1+round((STRK.w_z((j*(k-1))+k,1)/min(STRK.w_z))*99);
%             end
%             plot3(LWT.C((j*(k-1))+k,1),LWT.C((j*(k-1))+k,2),LWT.C((j*(k-1))+k,3),'o','Color',C_map(color_code_wz,:))
%             plot3(RWT.C((j*(k-1))+k,1),RWT.C((j*(k-1))+k,2),RWT.C((j*(k-1))+k,3),'o','Color',C_map(color_code_wz,:))
% 
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     colorbar
%     hold off
    
    
    
    
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             plot3(LWT.C((j*(k-1))+k,1),LWT.C((j*(k-1))+k,2),LWT.C((j*(k-1))+k,3),'o','Color','b')
%             plot3(RWT.C((j*(k-1))+k,1),RWT.C((j*(k-1))+k,2),RWT.C((j*(k-1))+k,3),'o','Color','b')
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             plot3(LWT.C((j*(k-1))+k,1),LWT.C((j*(k-1))+k,2),LWT.C((j*(k-1))+k,3),'o','Color','g')
%             plot3(RWT.C((j*(k-1))+k,1),RWT.C((j*(k-1))+k,2),RWT.C((j*(k-1))+k,3),'o','Color','g')
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold off
%     
%     figure()
%     hold on
%     for j = 1:nr_of_seq
%         for k = 1:nr_wb
%             plot3(LWT.D((j*(k-1))+k,1),LWT.D((j*(k-1))+k,2),LWT.D((j*(k-1))+k,3),'o','Color','m')
%             plot3(RWT.D((j*(k-1))+k,1),RWT.D((j*(k-1))+k,2),RWT.D((j*(k-1))+k,3),'o','Color','m')
%         end
%     end
%     surf(x,y,z,'FaceColor','black','EdgeColor','none');
%     alpha(0.2)
%     axis equal
%     hold off

end

