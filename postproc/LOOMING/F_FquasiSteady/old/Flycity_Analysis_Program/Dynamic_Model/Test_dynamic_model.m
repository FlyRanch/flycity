function Test_dynamic_model( settings, pathDB )


    % Program to test the dynamic model by generating a test-case as input
    % for the dynamical model:
    
    % Import an average wingbeat:
    
    seq_nr = 1;
    
    wb_id = 22;
    
    n_pol_theta = settings.n_pol_theta;
    n_pol_eta   = settings.n_pol_eta;
    n_pol_phi   = settings.n_pol_phi;
    
    a_theta_L = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_eta_L   = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_phi_L   = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
    
    a_theta_R = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_eta_R   = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_phi_R   = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
   
    down_up = pathDB.poly_fit.a_avg.down_up(seq_nr);
    
    f       = pathDB.poly_fit.a_avg.f(seq_nr);

%     a_theta_L   = pathDB.poly_fit.a_fit.theta_L(:,wb_id,seq_nr);
%     a_eta_L     = pathDB.poly_fit.a_fit.eta_L(:,wb_id,seq_nr);
%     a_phi_L     = pathDB.poly_fit.a_fit.phi_L(:,wb_id,seq_nr);
%     
%     a_theta_R   = pathDB.poly_fit.a_fit.theta_R(:,wb_id,seq_nr);
%     a_eta_R     = pathDB.poly_fit.a_fit.eta_R(:,wb_id,seq_nr);
%     a_phi_R     = pathDB.poly_fit.a_fit.phi_R(:,wb_id,seq_nr);
%     
%     down_up     = pathDB.poly_fit.a_fit.down_up(wb_id,seq_nr);
%     
%     f           = pathDB.poly_fit.a_fit.f(wb_id,seq_nr);
%     
%     wb_loc_filt = pathDB.wingbeats.wingbeat_loc(wb_id,:,seq_nr);
%     
%     theta_L_filt = pathDB.wing_kin.theta_L(seq_nr,wb_loc_filt(1):wb_loc_filt(2));
%     eta_L_filt   = pathDB.wing_kin.eta_L(seq_nr,wb_loc_filt(1):wb_loc_filt(2));
%     phi_L_filt   = pathDB.wing_kin.phi_L(seq_nr,wb_loc_filt(1):wb_loc_filt(2));
%     
%     theta_R_filt = pathDB.wing_kin.theta_R(seq_nr,wb_loc_filt(1):wb_loc_filt(2));
%     eta_R_filt   = pathDB.wing_kin.eta_R(seq_nr,wb_loc_filt(1):wb_loc_filt(2));
%     phi_R_filt   = pathDB.wing_kin.phi_R(seq_nr,wb_loc_filt(1):wb_loc_filt(2));
%     
%     wL_filt     = pathDB.filt.wL(wb_loc_filt(1):wb_loc_filt(2),:,seq_nr)';
%     wR_filt     = pathDB.filt.wR(wb_loc_filt(1):wb_loc_filt(2),:,seq_nr)';
    
    nr_points = 500;
    
%     nr_points = pathDB.wingbeats.wingbeat_duration(wb_id,seq_nr)+1;
        
    dt = (1/f)/(nr_points-1);
    
    a_fit_L.a_theta     = a_theta_L;
    a_fit_L.a_eta       = a_eta_L;
    a_fit_L.a_phi       = a_phi_L;
    a_fit_L.f           = f;
    a_fit_L.down_up     = down_up;
    
    a_fit_R.a_theta     = a_theta_R;
    a_fit_R.a_eta       = a_eta_R;
    a_fit_R.a_phi       = a_phi_R;
    a_fit_R.f           = f;
    a_fit_R.down_up     = down_up;
    
    [ kine_L ] = poly_derivatives( a_fit_L, nr_points );
    
    [ kine_R ] = poly_derivatives( a_fit_R, nr_points );
    
    t               = kine_L.t;
    theta_L         = kine_L.theta;
    eta_L           = kine_L.eta;
    phi_L           = kine_L.phi;
    theta_dot_L     = kine_L.theta_dot;
    eta_dot_L       = -kine_L.eta_dot;
    phi_dot_L       = -kine_L.phi_dot;
    theta_ddot_L    = kine_L.theta_ddot;
    eta_ddot_L      = -kine_L.eta_ddot;
    phi_ddot_L      = -kine_L.phi_ddot;
    theta_R         = kine_R.theta;
    eta_R           = kine_R.eta;
    phi_R           = kine_R.phi;
    theta_dot_R     = -kine_R.theta_dot;
    eta_dot_R       = -kine_R.eta_dot;
    phi_dot_R       = kine_R.phi_dot;
    theta_ddot_R    = -kine_R.theta_ddot;
    eta_ddot_R      = -kine_R.eta_ddot;
    phi_ddot_R      = kine_R.phi_ddot;
    
        
    figure()
    hold on
    subplot(3,1,1); plot(t,theta_L,'g',t,eta_L,'r',t,phi_L,'b')
    subplot(3,1,2); plot(t,theta_dot_L,'g',t,eta_dot_L,'r',t,phi_dot_L,'b')
    subplot(3,1,3); plot(t,theta_ddot_L,'g',t,eta_ddot_L,'r',t,phi_ddot_L,'b')
    hold off
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_L,'r',t(1:(end-1)),theta_L_filt,'k',t,theta_R,'g',t(1:(end-1)),theta_R_filt,'k')
%     subplot(3,1,2); plot(t,eta_L,'r',t(1:(end-1)),eta_L_filt,'k',t,eta_R,'g',t(1:(end-1)),eta_R_filt,'k')
%     subplot(3,1,3); plot(t,phi_L,'r',t(1:(end-1)),phi_L_filt,'k',t,phi_R,'g',t(1:(end-1)),phi_R_filt,'k')
%     hold off
        
    
    wb_loc      = [ 1 nr_points];
    down_loc    = [ 1 round(down_up*nr_points)];
    up_loc      = [ round(down_up*nr_points)+1 nr_points];
    
    R_strk = pathDB.rot_mat.Rstr;
    
    wing_kin.R_strk  = R_strk;
    wing_kin.theta_L = theta_L;
    wing_kin.eta_L   = eta_L;
    wing_kin.phi_L   = phi_L;
    wing_kin.theta_R = theta_R;
    wing_kin.eta_R   = eta_R;
    wing_kin.phi_R   = phi_R;
    wing_kin.theta_dot_L = theta_dot_L;
    wing_kin.eta_dot_L   = eta_dot_L;
    wing_kin.phi_dot_L   = phi_dot_L;
    wing_kin.theta_dot_R = theta_dot_R;
    wing_kin.eta_dot_R   = eta_dot_R;
    wing_kin.phi_dot_R   = phi_dot_R;
    wing_kin.theta_ddot_L = theta_ddot_L;
    wing_kin.eta_ddot_L   = eta_ddot_L;
    wing_kin.phi_ddot_L   = phi_ddot_L;
    wing_kin.theta_ddot_R = theta_ddot_R;
    wing_kin.eta_ddot_R   = eta_ddot_R;
    wing_kin.phi_ddot_R   = phi_ddot_R;
    
    [ rot_mat ] = wingkin_rotation_mat( wing_kin, nr_points );
    
    R_theta_L   = rot_mat.R_theta_L;
    R_eta_L     = rot_mat.R_eta_L;
    R_phi_L     = rot_mat.R_phi_L;
    R_theta_R   = rot_mat.R_theta_R;
    R_eta_R     = rot_mat.R_eta_R;
    R_phi_R     = rot_mat.R_phi_R;
    RL          = rot_mat.RL;
    RR          = rot_mat.RR;
    wL          = rot_mat.wL;
    wR          = rot_mat.wR;
    w_dot_L     = rot_mat.w_dot_L;
    w_dot_R     = rot_mat.w_dot_R;
    

    
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL(1,:),'r',t(1:(end-1)),wL_filt(1,:),'k')
%     subplot(3,1,2); plot(t,wL(2,:),'r',t(1:(end-1)),wL_filt(2,:),'k')
%     subplot(3,1,3); plot(t,wL(3,:),'r',t(1:(end-1)),wL_filt(3,:),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wR(1,:),'r',t(1:(end-1)),wR_filt(1,:),'k')
%     subplot(3,1,2); plot(t,wR(2,:),'r',t(1:(end-1)),wR_filt(2,:),'k')
%     subplot(3,1,3); plot(t,wR(3,:),'r',t(1:(end-1)),wR_filt(3,:),'k')
%     hold off
    
    kine.u_strk                 = zeros(3,nr_points);
    kine.w_strk                 = zeros(3,nr_points);
    kine.wL                     = wL;
    kine.wR                     = wR;
    kine.RL                     = RL;
    kine.RR                     = RR;
    kine.R_strk                 = R_strk;
    
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    
    wb.wb_loc                   = wb_loc;
    wb.down_loc                 = down_loc;
    wb.up_loc                   = up_loc;
    wb.dt                       = dt;
    
    [ FM_strkpln1, FM_L1, FM_R1 ,~, ~, ~, ~] = Aerodynamic_forces( kine, body_model, wing_model, wb, 0);
    
    [ FM_strkpln2, FM_L2, FM_R2 ,UL, UR, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R] = Aerodynamic_forces( kine, body_model, wing_model, wb, 1);
    
    figure()
    hold on
    subplot(2,1,1); plot(radtodeg(alfa_L(end,:)))
    subplot(2,1,2); plot(radtodeg(alfa_R(end,:)))
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(diff(alfa_L(end,:)))
    subplot(2,1,2); plot(diff(alfa_R(end,:)))
    hold off
    
    figure()
    hold on
    subplot(2,1,1); plot(1:nr_points,alfa_dot_L(end,:),'b',down_loc(2),0,'o')
    subplot(2,1,2); plot(1:nr_points,alfa_dot_R(end,:),'b',down_loc(2),0,'o')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,UL(1,:,end),'r',t,UR(1,:,end),'g')
    subplot(3,1,2); plot(t,UL(2,:,end),'r',t,UR(2,:,end),'g')
    subplot(3,1,3); plot(t,UL(3,:,end),'r',t,UR(3,:,end),'g')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,FM_strkpln1(1,:),'k',t,FM_strkpln2(1,:),'b')
    subplot(3,1,2); plot(t,FM_strkpln1(2,:),'k',t,FM_strkpln2(2,:),'b')
    subplot(3,1,3); plot(t,FM_strkpln1(3,:),'k',t,FM_strkpln2(3,:),'b')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,FM_strkpln1(4,:),'k',t,FM_strkpln2(4,:),'b')
    subplot(3,1,2); plot(t,FM_strkpln1(5,:),'k',t,FM_strkpln2(5,:),'b')
    subplot(3,1,3); plot(t,FM_strkpln1(6,:),'k',t,FM_strkpln2(6,:),'b')
    hold off
    
    mean(FM_strkpln1(3,:))
    
    mean(FM_strkpln2(3,:))
    
    pathDB.body_model.mass_fly(seq_nr)*9.81
    
    mean(FM_strkpln1(5,:))
    
    mean(FM_strkpln2(5,:))
    
    wtL1 = [0.05; -1; 0];
    wtR1 = [0.05; 1; 0];
    wtL2 = [-0.05; -1; 0];
    wtR2 = [-0.05; 1; 0];
    wtL3 = [0.025; -1; 0];
    wtR3 = [0.025; 1; 0];
    
    figure()
    hold on
    for i = 1:nr_points
        xyz1_L = RL(:,:,i)'*wtL1;
        xyz2_L = RL(:,:,i)'*wtL2;
        xyz1_R = RR(:,:,i)'*wtR1;
        xyz2_R = RR(:,:,i)'*wtR2;
        plot3([xyz1_L(1) xyz2_L(1)],[xyz1_L(2) xyz2_L(2)],[xyz1_L(3) xyz2_L(3)],'r')
        plot3([xyz1_R(1) xyz2_R(1)],[xyz1_R(2) xyz2_R(2)],[xyz1_R(3) xyz2_R(3)],'g')
    end
    axis equal
    hold off
    
    quiv_scale = 5000;
    
    figure()
    hold on
    for i = 1:nr_points
        xyz1_L  = R_strk*RL(:,:,i)'*wtL1;
        xyz2_L  = R_strk*RL(:,:,i)'*wtL2;
        xyz3_L  = R_strk*RL(:,:,i)'*wtL3;
        xyz1_R  = R_strk*RR(:,:,i)'*wtR1;
        xyz2_R  = R_strk*RR(:,:,i)'*wtR2;
        xyz3_R  = R_strk*RR(:,:,i)'*wtR3;       
        F_left  = R_strk*RL(:,:,i)'*(quiv_scale.*FM_L2(1:3,i));
        F_right = R_strk*RR(:,:,i)'*(quiv_scale.*FM_R2(1:3,i));
        plot3([xyz1_L(1) xyz2_L(1)],[xyz1_L(2) xyz2_L(2)],[xyz1_L(3) xyz2_L(3)],'k')
        plot3([xyz1_R(1) xyz2_R(1)],[xyz1_R(2) xyz2_R(2)],[xyz1_R(3) xyz2_R(3)],'k')
        plot3(xyz1_L(1),xyz1_L(2),xyz1_L(3),'.','Color','k')
        plot3(xyz1_R(1),xyz1_R(2),xyz1_R(3),'.','Color','k')
        quiver3(xyz3_L(1),xyz3_L(2),xyz3_L(3),F_left(1),F_left(2),F_left(3),'r')
        quiver3(xyz3_R(1),xyz3_R(2),xyz3_R(3),F_right(1),F_right(2),F_right(3),'g')
    end
    axis equal
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wL(1,:),'r',t,wR(1,:),'g')
    subplot(3,1,2); plot(t,wL(2,:),'r',t,wR(2,:),'g')
    subplot(3,1,3); plot(t,wL(3,:),'r',t,wR(3,:),'g')
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,w_dot_L(1,:),'r',t,w_dot_R(1,:),'g')
    subplot(3,1,2); plot(t,w_dot_L(2,:),'r',t,w_dot_R(2,:),'g')
    subplot(3,1,3); plot(t,w_dot_L(3,:),'r',t,w_dot_R(3,:),'g')
    hold off
    
end

