function Force_plot3(settings,pathDB,seq_nr,wing_l,qL1,qL2,qL3,qL4,qR1,qR2,qR3,qR4,u_wing_L,v_wing_L,w_wing_L,u_wing_R,v_wing_R,w_wing_R,alfa_wing_L,alfa_wing_R,Lwingbeat_loc,Rwingbeat_loc,Lwingtip,Rwingtip,left_on,right_on,t,fg_nr1, fg_nr2, fg_nr3,save_on_off)

    
    % Program that calculates the average generated force by the left and
    % the right wing an plots these forces for every upstroke/downstroke
    % and the combination of up and downstroke together with the average body
    % orientation and average wing orientation.
    
    % addpath body model
    
    addpath(char(settings.path_names(2)));
    addpath(char(settings.path_names(7)));
    addpath(char(settings.path_names(8)));
    
    % Load the wing model into the function and calculate the location of
    % the sections and the average chord of each section:
    
    N = length(u_wing_L(1,:));
    
    pos = [1 1-(0.5/(N-2)):-(1/(N-2)):(0.5/(N-2)) 0];
    
    x_sect = zeros(N-2,1);
    
    y_sect_L = -wing_l.*pos(2:N-1);
    
    y_sect_R = wing_l.*pos(2:N-1);
    
    delta_R = wing_l/(N-2);
    
    joint_pos_L = pathDB.joint_pos_L(:,seq_nr);
    
    joint_pos_R = pathDB.joint_pos_R(:,seq_nr);
    
    % Load wing and body model --------------------------------------------
    
    cd(char(settings.sequence_names(seq_nr)))
    
    cd('flytracks')
    
    load('ManualFit_flytracks');
    
    cd ..
    
    cd ..
    
    % Assign model parameters
    PAR.params = ManualFit.params;
    PAR.DLT = ManualFit.DLT;
    PAR.cam = ManualFit.cam;
    
    wing_scale = ManualFit.params.wingscale;
    
    body_scale = ManualFit.params.bodyscale;
    
    %clear ManualFit

    % Define the Tracking Parameters
    PAR.pixpermm = 1;
    PAR.numfly = 1;
    %Number of parameters of the model (i.e. 8 control points)
    PAR.mdlpar = 15*ones(1,PAR.numfly);
    PAR.statedim = PAR.mdlpar;
    PAR.modelfun_H = @modcurvesplineP;
    PAR.etamax = 0;
    
    %spline order
    PAR.c = 4;
    PAR.L1 = 15; %# of steps for body along length
    PAR.L2 = 6; %# of steps for head along length
    PAR.L3 = 25; %# of steps for wing around the boundary
    PAR.T1 = 13; %# of theta steps for head and body
    PAR.T2 = 2; %# of steps towards center of wing   
    
    SOLN = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1];
    
    clear flymodQ
    [x,y,z] = flymodQ(SOLN,PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    
%     figure()
%     plot(y{2},x{2})
%     axis equal
%     
%     figure()
%     plot(x{3}-joint_pos_R(1),y{3}-joint_pos_R(2))
%     hold on
%     scatter(x_sect,y_sect_R,'g')
%     plot(x{2}-joint_pos_L(1),y{2}-joint_pos_L(2))
%     scatter(x_sect,y_sect_L,'r')
%     axis equal
%     hold off
    
    % Extract chord length at the section centers and assume it is equal to
    % the average chord length for the section. This assumption is not very
    % accurate, however it will be sufficient if the number of sections is
    % big enough:
    
    x_wingR = x{3};
    y_wingR = y{3};
    
    chords = zeros(N-2,1);
    M = length(x_wingR(:,1));
    
    for i = 1:N-2
    int_val1 = interp1(y_wingR(1:((M+1)/2),1)-joint_pos_R(2),x_wingR(1:((M+1)/2),1)-joint_pos_R(1),y_sect_R(i));
    int_val2 = interp1(y_wingR(((M+1)/2):M,1)-joint_pos_R(2),x_wingR(((M+1)/2):M,1)-joint_pos_R(1),y_sect_R(i));
    
    chords(i) = int_val1-int_val2;
    end
    
    % Now calculate for every wingbeat the quasi-steady forces generated
    % during the up- and down-stroke for the left and right wingbeat:
    
        
    if Lwingbeat_loc(1,1) > Lwingbeat_loc(1,2)
        start_L = 0;
        end_L = length(Lwingbeat_loc(:,1))-1;
    else
        start_L = 1;
        end_L = length(Lwingbeat_loc(:,1))-2;
    end
    
    if Rwingbeat_loc(1,1) > Rwingbeat_loc(1,2)
        start_R = 0;
        end_R = length(Rwingbeat_loc(:,1))-1;
    else
        start_R = 1;
        end_R = length(Rwingbeat_loc(:,1))-2;
    end
    
    
    % Left wingbeat:

    
    % Downstroke:
    
    Fx_L_down = zeros(end_L,1);
    Fy_L_down = zeros(end_L,1);
    Fz_L_down = zeros(end_L,1);
    
    Fn_L_down = nan(end_L,50);
    Ft_L_down = nan(end_L,50);
    
    Mn_L_down = nan(end_L,50);
    Mt_L_down = nan(end_L,50);
    
    R_F_L_down = nan(end_L,50);
    
    Mx_L_down = zeros(end_L,1);
    My_L_down = zeros(end_L,1);
    Mz_L_down = zeros(end_L,1);
    
    q_avg_L_down = zeros(end_L,4);
    
    r_body_L_down = zeros(3,end_L);
    
    % Upstroke:
    
    Fx_L_up = zeros(end_L,1);
    Fy_L_up = zeros(end_L,1);
    Fz_L_up = zeros(end_L,1);
    
    Fn_L_up = nan(end_L,50);
    Ft_L_up = nan(end_L,50);
    
    Mn_L_up = nan(end_L,50);
    Mt_L_up = nan(end_L,50);
    
    R_F_L_up = nan(end_L,50);
    
    Mx_L_up = zeros(end_L,1);
    My_L_up = zeros(end_L,1);
    Mz_L_up = zeros(end_L,1);
    
    q_avg_L_up = zeros(end_L,4);
    
    r_body_L_up = zeros(3,end_L);
    
    
    for i = 1:end_L
        
        y_sect = y_sect_L;
        
        % Downstroke:
        
        a = Lwingbeat_loc(i,2):1:(Lwingbeat_loc(i+start_L,1)-1);
        
        q1 = qL1(a);
        q2 = qL2(a);
        q3 = qL3(a);
        q4 = qL4(a);
        
        u_wing = u_wing_L(a,:);
        v_wing = v_wing_L(a,:);
        w_wing = w_wing_L(a,:);
        
        alfa_wing = alfa_wing_L(a,:);
        
        up = 0;
        down = 1;

        [Fx_L_down(i),Fy_L_down(i),Fz_L_down(i),Fn_L_down(i,:),Ft_L_down(i,:),Mn_L_down(i,:),Mt_L_down(i,:),R_F_L_down(i,:),r_body_L_down(:,i),Mx_L_down(i,1),My_L_down(i,1),Mz_L_down(i,1),q_avg_L_down(i,:)] = Quasi_steady_force(q1,q2,q3,q4,u_wing,v_wing,w_wing,alfa_wing,chords,delta_R,y_sect,up,down);
        
        clear a b up down q1 q2 q3 q4 u_wing v_wing w_wing alfa_wing
        
        % Upstroke:
        
        b = Lwingbeat_loc(i+start_L,1):1:(Lwingbeat_loc(i+1,2)-1);
        
        q1 = qL1(b);
        q2 = qL2(b);
        q3 = qL3(b);
        q4 = qL4(b);
        
        u_wing = u_wing_L(b,:);
        v_wing = v_wing_L(b,:);
        w_wing = w_wing_L(b,:);
        
        alfa_wing = alfa_wing_L(b,:);
        
        up = 1;
        down = 0;
        
        [Fx_L_up(i),Fy_L_up(i),Fz_L_up(i),Fn_L_up(i,:),Ft_L_up(i,:),Mn_L_up(i,:),Mt_L_up(i,:),R_F_L_up(i,:),r_body_L_up(:,i),Mx_L_up(i,1),My_L_up(i,1),Mz_L_up(i,1),q_avg_L_up(i,:)] = Quasi_steady_force(q1,q2,q3,q4,u_wing,v_wing,w_wing,alfa_wing,chords,delta_R,y_sect,up,down);
       
        
        clear a b up down q1 q2 q3 q4 u_wing v_wing w_wing alfa_wing
    end
      
    
    % Right wingbeat:

    
    % Downstroke:
    
    Fx_R_down = zeros(end_R,1);
    Fy_R_down = zeros(end_R,1);
    Fz_R_down = zeros(end_R,1);
    
    Fn_R_down = nan(end_R,50);
    Ft_R_down = nan(end_R,50);
    
    Mn_R_down = nan(end_R,50);
    Mt_R_down = nan(end_R,50);
    
    R_F_R_down = nan(end_R,50);
    
    Mx_R_down = zeros(end_R,1);
    My_R_down = zeros(end_R,1);
    Mz_R_down = zeros(end_R,1);
    
    q_avg_R_down = zeros(end_R,4);
    
    r_body_R_down = zeros(3,end_R);
    
    % Upstroke:
    
    Fx_R_up = zeros(end_R,1);
    Fy_R_up = zeros(end_R,1);
    Fz_R_up = zeros(end_R,1);
    
    Fn_R_up = nan(end_R,50);
    Ft_R_up = nan(end_R,50);
    
    Mn_R_up = nan(end_R,50);
    Mt_R_up = nan(end_R,50);
    
    R_F_R_up = nan(end_R,50);
    
    Mx_R_up = zeros(end_R,1);
    My_R_up = zeros(end_R,1);
    Mz_R_up = zeros(end_R,1);
    
    q_avg_R_up = zeros(end_R,4);
    
    r_body_R_up = zeros(3,end_R);
    
    for i = 1:end_R
        
        y_sect = y_sect_R;
        
        % Downstroke:
        
        a = Rwingbeat_loc(i,2):1:(Rwingbeat_loc(i+start_R,1)-1);
        
        q1 = qR1(a);
        q2 = qR2(a);
        q3 = qR3(a);
        q4 = qR4(a);
        
        u_wing = u_wing_R(a,:);
        v_wing = v_wing_R(a,:);
        w_wing = w_wing_R(a,:);
        
        alfa_wing = alfa_wing_R(a,:);
        
        up = 0;
        down = 1;

        [Fx_R_down(i),Fy_R_down(i),Fz_R_down(i),Fn_R_down(i,:),Ft_R_down(i,:),Mn_R_down(i,:),Mt_R_down(i,:),R_F_R_down(i,:),r_body_R_down(:,i),Mx_R_down(i,1),My_R_down(i,1),Mz_R_down(i,1),q_avg_R_down(i,:)] = Quasi_steady_force(q1,q2,q3,q4,u_wing,v_wing,w_wing,alfa_wing,chords,delta_R,y_sect,up,down);
        
        clear a b up down q1 q2 q3 q4 u_wing v_wing w_wing alfa_wing
        
        % Upstroke:
        
        b = Rwingbeat_loc(i+start_R,1):1:(Rwingbeat_loc(i+1,2)-1);
        
        q1 = qR1(b);
        q2 = qR2(b);
        q3 = qR3(b);
        q4 = qR4(b);
        
        u_wing = u_wing_R(b,:);
        v_wing = v_wing_R(b,:);
        w_wing = w_wing_R(b,:);
        
        alfa_wing = alfa_wing_R(b,:);
        
        up = 1;
        down = 0;
        
        [Fx_R_up(i),Fy_R_up(i),Fz_R_up(i),Fn_R_up(i,:),Ft_R_up(i,:),Mn_R_up(i,:),Mt_R_up(i,:),R_F_R_up(i,:),r_body_R_up(:,i),Mx_R_up(i,1),My_R_up(i,1),Mz_R_up(i,1),q_avg_R_up(i,:)] = Quasi_steady_force(q1,q2,q3,q4,u_wing,v_wing,w_wing,alfa_wing,chords,delta_R,y_sect,up,down);
       
        
        clear a b up down q1 q2 q3 q4 u_wing v_wing w_wing alfa_wing
    end
    
    % Plot the development of the forces during an upstroke and a downstroke

    
    for i = 1:end_R
        
        
        % Downstroke ------------------------------------------------------
        
        a_L = Lwingbeat_loc(i,2):1:(Lwingbeat_loc(i+start_L,1)-1);
        
        q1 = qL1(a_L);
        q2 = qL2(a_L);
        q3 = qL3(a_L);
        q4 = qL4(a_L);
        
        N_a_L = length(a_L);
                
        L_wt_down = Lwingtip(a_L,:);
        
        LE_L_down = zeros(3,N_a_L);
        TE_L_down = zeros(3,N_a_L);
    
        LE_R_down = zeros(3,N_a_L);
        TE_R_down = zeros(3,N_a_L);
        
        
        for j = 1:N_a_L
                        
            DCM = quat2matNEW([q1(j) q2(j) q3(j) q4(j)]);
            
            LE_L_down(:,j) = pathDB.wing_l(seq_nr).*DCM*[0.03; -1; 0];
            TE_L_down(:,j) = pathDB.wing_l(seq_nr).*DCM*[-0.05; -1; 0];
            
            clear DCM 
        end
        
        a_R = Rwingbeat_loc(i,2):1:(Rwingbeat_loc(i+start_R,1)-1);
        
        q1 = qR1(a_R);
        q2 = qR2(a_R);
        q3 = qR3(a_R);
        q4 = qR4(a_R);
        
        N_a_R = length(a_R);
        
        R_wt_down = Rwingtip(a_R,:);
        
        for j = 1:N_a_R
            
            DCM = quat2matNEW([q1(j) q2(j) q3(j) q4(j)]);
                        
            LE_R_down(:,j) = pathDB.wing_l(seq_nr).*DCM*[0.03; 1; 0];
            TE_R_down(:,j) = pathDB.wing_l(seq_nr).*DCM*[-0.05; 1; 0];
            
            clear DCM
        end

        
        
        % Upstroke --------------------------------------------------------
        
        b_L = Lwingbeat_loc(i+start_L,1):1:(Lwingbeat_loc(i+1,2)-1);
        
        q1 = qL1(b_L);
        q2 = qL2(b_L);
        q3 = qL3(b_L);
        q4 = qL4(b_L);
        
        N_b_L = length(b_L);

        
        L_wt_up = Lwingtip(b_L,:);

        
        LE_L_up = zeros(3,N_b_L);
        TE_L_up = zeros(3,N_b_L);
    
        LE_R_up = zeros(3,N_b_L);
        TE_R_up = zeros(3,N_b_L);
        
        for j = 1:N_b_L
            
            
            DCM = quat2matNEW([q1(j) q2(j) q3(j) q4(j)]);
            
            
            LE_L_up(:,j) = pathDB.wing_l(seq_nr).*DCM*[0.03; -1; 0];
            TE_L_up(:,j) = pathDB.wing_l(seq_nr).*DCM*[-0.05; -1; 0];
            
            clear DCM 
        end
        
        b_R = Rwingbeat_loc(i+start_R,1):1:(Rwingbeat_loc(i+1,2)-1);
        
        q1 = qR1(b_R);
        q2 = qR2(b_R);
        q3 = qR3(b_R);
        q4 = qR4(b_R);
        
        N_b_R = length(b_R);

        
        R_wt_up = Rwingtip(b_R,:);

        
        for j = 1:N_b_R
            

            
            DCM = quat2matNEW([q1(j) q2(j) q3(j) q4(j)]);
            
            
            LE_R_up(:,j) = pathDB.wing_l(seq_nr).*DCM*[0.03; 1; 0];
            TE_R_up(:,j) = pathDB.wing_l(seq_nr).*DCM*[-0.05; 1; 0];
            
            clear DCM 
        end

        figure()
        surf(x{1},y{1},z{1},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong');
        hold on
        if left_on == 1
        %plot3(joint_pos_L(1)+L_wt_down(:,1),joint_pos_L(2)+L_wt_down(:,2),joint_pos_L(3)+L_wt_down(:,3),'k')
        %for k = 1:N_a_L
        %plot3([joint_pos_L(1)+LE_L_down(1,k) joint_pos_L(1)+TE_L_down(1,k)],[joint_pos_L(2)+LE_L_down(2,k) joint_pos_L(2)+TE_L_down(2,k)],[joint_pos_L(3)+LE_L_down(3,k) joint_pos_L(3)+TE_L_down(3,k)],'r')
        %end
        %plot3([joint_pos_L(1) joint_pos_L(1)+L_wt_down(1,1)],[joint_pos_L(2) joint_pos_L(2)+L_wt_down(1,2)],[joint_pos_L(3) joint_pos_L(3)+L_wt_down(1,3)],'k')
        %plot3([joint_pos_L(1) joint_pos_L(1)+L_wt_down(N_a_L,1)],[joint_pos_L(2) joint_pos_L(2)+L_wt_down(N_a_L,2)],[joint_pos_L(3) joint_pos_L(3)+L_wt_down(N_a_L,3)],'k')
        plot3([joint_pos_L(1) joint_pos_L(1)+r_body_L_down(1,i)],[joint_pos_L(2) joint_pos_L(2)+r_body_L_down(2,i)],[joint_pos_L(3) joint_pos_L(3)+r_body_L_down(3,i)],'k')
        quiver3(joint_pos_L(1)+r_body_L_down(1,i),joint_pos_L(2)+r_body_L_down(2,i),joint_pos_L(3)+r_body_L_down(3,i),1e5*Fx_L_down(i),1e5*Fy_L_down(i),1e5*Fz_L_down(i),'r')
        hold on
        end
        if right_on == 1
        %plot3(joint_pos_R(1)+R_wt_down(:,1),joint_pos_R(2)+R_wt_down(:,2),joint_pos_R(3)+R_wt_down(:,3),'k')
        %for k = 1:N_a_R
        %plot3([joint_pos_R(1)+LE_R_down(1,k) joint_pos_R(1)+TE_R_down(1,k)],[joint_pos_R(2)+LE_R_down(2,k) joint_pos_R(2)+TE_R_down(2,k)],[joint_pos_R(3)+LE_R_down(3,k) joint_pos_R(3)+TE_R_down(3,k)],'r')
        %end
        %plot3([joint_pos_R(1) joint_pos_R(1)+R_wt_down(1,1)],[joint_pos_R(2) joint_pos_R(2)+R_wt_down(1,2)],[joint_pos_R(3) joint_pos_R(3)+R_wt_down(1,3)],'k')
        %plot3([joint_pos_R(1) joint_pos_R(1)+R_wt_down(N_a_R,1)],[joint_pos_R(2) joint_pos_R(2)+R_wt_down(N_a_R,2)],[joint_pos_R(3) joint_pos_R(3)+R_wt_down(N_a_R,3)],'k')
        plot3([joint_pos_R(1) joint_pos_R(1)+r_body_R_down(1,i)],[joint_pos_R(2) joint_pos_R(2)+r_body_R_down(2,i)],[joint_pos_R(3) joint_pos_R(3)+r_body_R_down(3,i)],'k')
        quiver3(joint_pos_R(1)+r_body_R_down(1,i),joint_pos_R(2)+r_body_R_down(2,i),joint_pos_R(3)+r_body_R_down(3,i),1e5*Fx_R_down(i),1e5*Fy_R_down(i),1e5*Fz_R_down(i),'r')
        end
        hold on
        if left_on == 1
        hold on
        %plot3(joint_pos_L(1)+L_wt_up(:,1),joint_pos_L(2)+L_wt_up(:,2),joint_pos_L(3)+L_wt_up(:,3),'k')
        %for k = 1:N_b_L
        %plot3([joint_pos_L(1)+LE_L_up(1,k) joint_pos_L(1)+TE_L_up(1,k)],[joint_pos_L(2)+LE_L_up(2,k) joint_pos_L(2)+TE_L_up(2,k)],[joint_pos_L(3)+LE_L_up(3,k) joint_pos_L(3)+TE_L_up(3,k)],'g')
        %end
        %plot3([joint_pos_L(1) joint_pos_L(1)+L_wt_up(1,1)],[joint_pos_L(2) joint_pos_L(2)+L_wt_up(1,2)],[joint_pos_L(3) joint_pos_L(3)+L_wt_up(1,3)],'k')
        %plot3([joint_pos_L(1) joint_pos_L(1)+L_wt_up(N_b_L,1)],[joint_pos_L(2) joint_pos_L(2)+L_wt_up(N_b_L,2)],[joint_pos_L(3) joint_pos_L(3)+L_wt_up(N_b_L,3)],'k')
        plot3([joint_pos_L(1) joint_pos_L(1)+r_body_L_up(1,i)],[joint_pos_L(2) joint_pos_L(2)+r_body_L_up(2,i)],[joint_pos_L(3) joint_pos_L(3)+r_body_L_up(3,i)],'k')
        quiver3(joint_pos_L(1)+r_body_L_up(1,i),joint_pos_L(2)+r_body_L_up(2,i),joint_pos_L(3)+r_body_L_up(3,i),1e5*Fx_L_up(i),1e5*Fy_L_up(i),1e5*Fz_L_up(i),'g')
        hold on
        end
        if right_on == 1
        %plot3(joint_pos_R(1)+R_wt_up(:,1),joint_pos_R(2)+R_wt_up(:,2),joint_pos_R(3)+R_wt_up(:,3),'k')
        %for k = 1:N_b_R
        %plot3([joint_pos_R(1)+LE_R_up(1,k) joint_pos_R(1)+TE_R_up(1,k)],[joint_pos_R(2)+LE_R_up(2,k) joint_pos_R(2)+TE_R_up(2,k)],[joint_pos_R(3)+LE_R_up(3,k) joint_pos_R(3)+TE_R_up(3,k)],'g')
        %end
        %plot3([joint_pos_R(1) joint_pos_R(1)+R_wt_up(1,1)],[joint_pos_R(2) joint_pos_R(2)+R_wt_up(1,2)],[joint_pos_R(3) joint_pos_R(3)+R_wt_up(1,3)],'k')
        %plot3([joint_pos_R(1) joint_pos_R(1)+R_wt_up(N_b_R,1)],[joint_pos_R(2) joint_pos_R(2)+R_wt_up(N_b_R,2)],[joint_pos_R(3) joint_pos_R(3)+R_wt_up(N_b_R,3)],'k')
        plot3([joint_pos_R(1) joint_pos_R(1)+r_body_R_up(1,i)],[joint_pos_R(2) joint_pos_R(2)+r_body_R_up(2,i)],[joint_pos_R(3) joint_pos_R(3)+r_body_R_up(3,i)],'k')
        quiver3(joint_pos_R(1)+r_body_R_up(1,i),joint_pos_R(2)+r_body_R_up(2,i),joint_pos_R(3)+r_body_R_up(3,i),1e5*Fx_R_up(i),1e5*Fy_R_up(i),1e5*Fz_R_up(i),'g')
        end
        axis equal
        hold off
        
        
        

    end

 
   

    
    %----------------------------------------------------------------------
end
