function Conservation_test(settings,pathDB)

    % Test conservation of linear and angular momentum:
    
    seq_nr      = 1;
    nr_points   = 50;
    nr_wb       = 5;
    
    % Setup the body_model and wing_model:
    
    body_model.mass_fly         = pathDB.body_model.mass_fly(seq_nr);
    body_model.mass_body        = pathDB.body_model.mass_body(seq_nr);
    body_model.Joint_left       = pathDB.body_model.Joint_left(seq_nr,:)';
    body_model.Joint_right      = pathDB.body_model.Joint_right(seq_nr,:)';
    body_model.cg_b             = pathDB.body_model.cg(seq_nr,:)';
    body_model.Inertia          = pathDB.body_model.Inertia(:,:,seq_nr);
    body_model.g                = 1e-3*settings.g;
    
    wing_model.virtual_mass     = pathDB.wing_model.virtual_mass(seq_nr);
    wing_model.wing_cg_L        = pathDB.wing_model.wing_cg_L(seq_nr,:)';
    wing_model.wing_cg_R        = pathDB.wing_model.wing_cg_R(seq_nr,:)';
    wing_model.virtual_Inertia  = pathDB.wing_model.virtual_Inertia(:,:,seq_nr);
    wing_model.y_sect_L         = pathDB.wing_model.y_sect_L(:,:,seq_nr)';
    wing_model.chords_L         = pathDB.wing_model.chords_L(seq_nr,:)';
    wing_model.y_sect_R         = pathDB.wing_model.y_sect_R(:,:,seq_nr)';
    wing_model.chords_R         = pathDB.wing_model.chords_R(seq_nr,:)';
    wing_model.rho              = settings.rho_air;
    
    nr_sect                     = length(wing_model.y_sect_L(1,:));
    
    a_avg_theta_L   = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_avg_eta_L     = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_avg_phi_L     = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
    
    a_avg_theta_R   = pathDB.poly_fit.a_avg.theta_LR(:,seq_nr);
    a_avg_eta_R     = pathDB.poly_fit.a_avg.eta_LR(:,seq_nr);
    a_avg_phi_R     = pathDB.poly_fit.a_avg.phi_LR(:,seq_nr);  
   
    down_up         = pathDB.poly_fit.a_avg.down_up(seq_nr);
    
    f               = pathDB.poly_fit.a_avg.f(seq_nr);
    
    dt              = (1/f)/(nr_points-1);
    
    R_strk = pathDB.rot_mat.Rstr;
    
    n_pol_theta     = (length(a_avg_theta_L)-2)/2;
    n_pol_eta       = (length(a_avg_eta_L)-2)/2;
    n_pol_phi       = (length(a_avg_phi_L)-2)/2;
    
    % Set up the polynomial coefficients per wingbeat:
    
    a_theta_L       = nan(2*(n_pol_theta+1),nr_wb);
    a_eta_L         = nan(2*(n_pol_eta+1),nr_wb);
    a_phi_L         = nan(2*(n_pol_phi+1),nr_wb);
    a_theta_R       = nan(2*(n_pol_theta+1),nr_wb);
    a_eta_R         = nan(2*(n_pol_eta+1),nr_wb);
    a_phi_R         = nan(2*(n_pol_phi+1),nr_wb);
    
    for i = 1:nr_wb
        
        a_theta_L(:,i)  = a_avg_theta_L;
        a_eta_L(:,i)    = a_avg_eta_L;
        a_phi_L(:,i)    = a_avg_phi_L;
        a_theta_R(:,i)  = a_avg_theta_R;
        a_eta_R(:,i)    = a_avg_eta_R;
        a_phi_R(:,i)    = a_avg_phi_R;
        
    end
    
        
    % Create the wing kinematics:
    
    t               = nan(1,(nr_points*2)*nr_wb+1);
    theta_L         = nan(1,(nr_points*2)*nr_wb+1);
    eta_L           = nan(1,(nr_points*2)*nr_wb+1);
    phi_L           = nan(1,(nr_points*2)*nr_wb+1);
    theta_dot_L     = nan(1,(nr_points*2)*nr_wb+1);
    eta_dot_L       = nan(1,(nr_points*2)*nr_wb+1);
    phi_dot_L       = nan(1,(nr_points*2)*nr_wb+1);
    theta_ddot_L    = nan(1,(nr_points*2)*nr_wb+1);
    eta_ddot_L      = nan(1,(nr_points*2)*nr_wb+1);
    phi_ddot_L      = nan(1,(nr_points*2)*nr_wb+1);
    theta_R         = nan(1,(nr_points*2)*nr_wb+1);
    eta_R           = nan(1,(nr_points*2)*nr_wb+1);
    phi_R           = nan(1,(nr_points*2)*nr_wb+1);
    theta_dot_R     = nan(1,(nr_points*2)*nr_wb+1);
    eta_dot_R       = nan(1,(nr_points*2)*nr_wb+1);
    phi_dot_R       = nan(1,(nr_points*2)*nr_wb+1);
    theta_ddot_R    = nan(1,(nr_points*2)*nr_wb+1);
    eta_ddot_R      = nan(1,(nr_points*2)*nr_wb+1);
    phi_ddot_R      = nan(1,(nr_points*2)*nr_wb+1);
    RL              = nan(3,3,(nr_points*2)*nr_wb+1);
    RR              = nan(3,3,(nr_points*2)*nr_wb+1);
    wL              = nan(3,(nr_points*2)*nr_wb+1);
    wR              = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_L         = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_R         = nan(3,(nr_points*2)*nr_wb+1);
    wL_b            = nan(3,(nr_points*2)*nr_wb+1);
    wR_b            = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_L_b       = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_R_b       = nan(3,(nr_points*2)*nr_wb+1);
    wL_strk         = nan(3,(nr_points*2)*nr_wb+1);
    wR_strk         = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_L_strk    = nan(3,(nr_points*2)*nr_wb+1);
    w_dot_R_strk    = nan(3,(nr_points*2)*nr_wb+1);
    
    wb_loc          = nan(nr_wb,2);
    down_loc        = nan(nr_wb,2);
    up_loc          = nan(nr_wb,2);
    
        
    for i = 1:nr_wb
        
        a_fit.a_theta_L     = a_theta_L(:,i);
        a_fit.a_eta_L       = a_eta_L(:,i);
        a_fit.a_phi_L       = a_phi_L(:,i);
        a_fit.a_theta_R     = a_theta_R(:,i);
        a_fit.a_eta_R       = a_eta_R(:,i);
        a_fit.a_phi_R       = a_phi_R(:,i);
        a_fit.f             = f;
        a_fit.down_up       = down_up;
        a_fit.nr_points     = nr_points*2+1;
        a_fit.R_strk        = R_strk;

        [ kine ] = angular_velocities( a_fit );
        
        if i == 1
            
            range_i         = ((i-1)*(2*nr_points)+1):(i*(2*nr_points)+1);
            range_temp      = 1:(2*nr_points+1);
            
            t(range_i)              = kine.t(range_temp);
            theta_L(range_i)        = kine.theta_L(range_temp);
            eta_L(range_i)          = kine.eta_L(range_temp);
            phi_L(range_i)          = kine.phi_L(range_temp);
            theta_dot_L(range_i)    = kine.theta_dot_L(range_temp);
            eta_dot_L(range_i)      = kine.eta_dot_L(range_temp);
            phi_dot_L(range_i)      = kine.phi_dot_L(range_temp);
            theta_ddot_L(range_i)   = kine.theta_ddot_L(range_temp);
            eta_ddot_L(range_i)     = kine.eta_ddot_L(range_temp);
            phi_ddot_L(range_i)     = kine.phi_ddot_L(range_temp);
            theta_R(range_i)        = kine.theta_R(range_temp);
            eta_R(range_i)          = kine.eta_R(range_temp);
            phi_R(range_i)          = kine.phi_R(range_temp);
            theta_dot_R(range_i)    = kine.theta_dot_R(range_temp);
            eta_dot_R(range_i)      = kine.eta_dot_R(range_temp);
            phi_dot_R(range_i)      = kine.phi_dot_R(range_temp);
            theta_ddot_R(range_i)   = kine.theta_ddot_R(range_temp);
            eta_ddot_R(range_i)     = kine.eta_ddot_R(range_temp);
            phi_ddot_R(range_i)     = kine.phi_ddot_R(range_temp);
            RL(:,:,range_i)         = kine.RL(:,:,range_temp);
            RR(:,:,range_i)         = kine.RR(:,:,range_temp);
            wL(:,range_i)           = kine.wL(:,range_temp);
            wR(:,range_i)           = kine.wR(:,range_temp);
            w_dot_L(:,range_i)      = kine.w_dot_L(:,range_temp);
            w_dot_R(:,range_i)      = kine.w_dot_R(:,range_temp);
            wL_b(:,range_i)         = kine.wL_b(:,range_temp);
            wR_b(:,range_i)         = kine.wR_b(:,range_temp);
            w_dot_L_b(:,range_i)    = kine.w_dot_L_b(:,range_temp);
            w_dot_R_b(:,range_i)    = kine.w_dot_R_b(:,range_temp);
            wL_strk(:,range_i)      = kine.wL_strk(:,range_temp);
            wR_strk(:,range_i)      = kine.wR_strk(:,range_temp);
            w_dot_L_strk(:,range_i) = kine.w_dot_L_strk(:,range_temp);
            w_dot_R_strk(:,range_i) = kine.w_dot_R_strk(:,range_temp);
            
        elseif i == nr_wb && i > 1
            
            range_i         = ((i-1)*(2*nr_points)+1):(i*(2*nr_points)+1);
            range_temp      = 1:(2*nr_points+1);
            
            t(range_i)              = t((i-1)*(2*nr_points))+kine.dt+kine.t(range_temp);
            theta_L(range_i)        = kine.theta_L(range_temp);
            eta_L(range_i)          = kine.eta_L(range_temp);
            phi_L(range_i)          = kine.phi_L(range_temp);
            theta_dot_L(range_i)    = kine.theta_dot_L(range_temp);
            eta_dot_L(range_i)      = kine.eta_dot_L(range_temp);
            phi_dot_L(range_i)      = kine.phi_dot_L(range_temp);
            theta_ddot_L(range_i)   = kine.theta_ddot_L(range_temp);
            eta_ddot_L(range_i)     = kine.eta_ddot_L(range_temp);
            phi_ddot_L(range_i)     = kine.phi_ddot_L(range_temp);
            theta_R(range_i)        = kine.theta_R(range_temp);
            eta_R(range_i)          = kine.eta_R(range_temp);
            phi_R(range_i)          = kine.phi_R(range_temp);
            theta_dot_R(range_i)    = kine.theta_dot_R(range_temp);
            eta_dot_R(range_i)      = kine.eta_dot_R(range_temp);
            phi_dot_R(range_i)      = kine.phi_dot_R(range_temp);
            theta_ddot_R(range_i)   = kine.theta_ddot_R(range_temp);
            eta_ddot_R(range_i)     = kine.eta_ddot_R(range_temp);
            phi_ddot_R(range_i)     = kine.phi_ddot_R(range_temp);
            RL(:,:,range_i)         = kine.RL(:,:,range_temp);
            RR(:,:,range_i)         = kine.RR(:,:,range_temp);
            wL(:,range_i)           = kine.wL(:,range_temp);
            wR(:,range_i)           = kine.wR(:,range_temp);
            w_dot_L(:,range_i)      = kine.w_dot_L(:,range_temp);
            w_dot_R(:,range_i)      = kine.w_dot_R(:,range_temp);
            wL_b(:,range_i)         = kine.wL_b(:,range_temp);
            wR_b(:,range_i)         = kine.wR_b(:,range_temp);
            w_dot_L_b(:,range_i)    = kine.w_dot_L_b(:,range_temp);
            w_dot_R_b(:,range_i)    = kine.w_dot_R_b(:,range_temp);
            wL_strk(:,range_i)      = kine.wL_strk(:,range_temp);
            wR_strk(:,range_i)      = kine.wR_strk(:,range_temp);
            w_dot_L_strk(:,range_i) = kine.w_dot_L_strk(:,range_temp);
            w_dot_R_strk(:,range_i) = kine.w_dot_R_strk(:,range_temp);
           
        else
            
            range_i         = ((i-1)*(2*nr_points)+1):(i*(2*nr_points)+1);
            
            t(range_i)              = t((i-1)*(2*nr_points))+kine.dt+kine.t;
            theta_L(range_i)        = kine.theta_L;
            eta_L(range_i)          = kine.eta_L;
            phi_L(range_i)          = kine.phi_L;
            theta_dot_L(range_i)    = kine.theta_dot_L;
            eta_dot_L(range_i)      = kine.eta_dot_L;
            phi_dot_L(range_i)      = kine.phi_dot_L;
            theta_ddot_L(range_i)   = kine.theta_ddot_L;
            eta_ddot_L(range_i)     = kine.eta_ddot_L;
            phi_ddot_L(range_i)     = kine.phi_ddot_L;
            theta_R(range_i)        = kine.theta_R;
            eta_R(range_i)          = kine.eta_R;
            phi_R(range_i)          = kine.phi_R;
            theta_dot_R(range_i)    = kine.theta_dot_R;
            eta_dot_R(range_i)      = kine.eta_dot_R;
            phi_dot_R(range_i)      = kine.phi_dot_R;
            theta_ddot_R(range_i)   = kine.theta_ddot_R;
            eta_ddot_R(range_i)     = kine.eta_ddot_R;
            phi_ddot_R(range_i)     = kine.phi_ddot_R;
            RL(:,:,range_i)         = kine.RL;
            RR(:,:,range_i)         = kine.RR;
            wL(:,range_i)           = kine.wL;
            wR(:,range_i)           = kine.wR;
            w_dot_L(:,range_i)      = kine.w_dot_L;
            w_dot_R(:,range_i)      = kine.w_dot_R;
            wL_b(:,range_i)         = kine.wL_b;
            wR_b(:,range_i)         = kine.wR_b;
            w_dot_L_b(:,range_i)    = kine.w_dot_L_b;
            w_dot_R_b(:,range_i)    = kine.w_dot_R_b;
            wL_strk(:,range_i)      = kine.wL_strk;
            wR_strk(:,range_i)      = kine.wR_strk;
            w_dot_L_strk(:,range_i) = kine.w_dot_L_strk;
            w_dot_R_strk(:,range_i) = kine.w_dot_R_strk;
            
        end
        
        wb_loc(i,:)         = [ ((i-1)*nr_points+1) (i*nr_points) ];
        down_loc(i,:)       = [ ((i-1)*nr_points+1) (((i-1)*nr_points)+round(down_up*nr_points)) ];
        up_loc(i,:)         = [ (((i-1)*nr_points)+round(down_up*nr_points)+1) (i*nr_points) ];
            
        clear kine a_fit
        
    end
    
    figure()
    hold on
    subplot(3,1,1); plot(t,theta_L,t,ones(((nr_points*2)*nr_wb+1),1)*mean(theta_L))
    subplot(3,1,2); plot(t,eta_L,t,ones(((nr_points*2)*nr_wb+1),1)*mean(eta_L))
    subplot(3,1,3); plot(t,phi_L,t,ones(((nr_points*2)*nr_wb+1),1)*mean(phi_L))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,theta_R,t,ones(((nr_points*2)*nr_wb+1),1)*mean(theta_R))
    subplot(3,1,2); plot(t,eta_R,t,ones(((nr_points*2)*nr_wb+1),1)*mean(eta_R))
    subplot(3,1,3); plot(t,phi_R,t,ones(((nr_points*2)*nr_wb+1),1)*mean(phi_R))
    hold off
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_L,t,ones(((nr_points*2)*nr_wb),1)*mean(theta_dot_L))
%     subplot(3,1,2); plot(t,eta_dot_L,t,ones(((nr_points*2)*nr_wb),1)*mean(eta_dot_L))
%     subplot(3,1,3); plot(t,phi_dot_L,t,ones(((nr_points*2)*nr_wb),1)*mean(phi_dot_L))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_R,t,ones(((nr_points*2)*nr_wb),1)*mean(theta_dot_R))
%     subplot(3,1,2); plot(t,eta_dot_R,t,ones(((nr_points*2)*nr_wb),1)*mean(eta_dot_R))
%     subplot(3,1,3); plot(t,phi_dot_R,t,ones(((nr_points*2)*nr_wb),1)*mean(phi_dot_R))
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_L,t,ones(((nr_points*2)*nr_wb),1)*mean(theta_ddot_L))
%     subplot(3,1,2); plot(t,eta_ddot_L,t,ones(((nr_points*2)*nr_wb),1)*mean(eta_ddot_L))
%     subplot(3,1,3); plot(t,phi_ddot_L,t,ones(((nr_points*2)*nr_wb),1)*mean(phi_ddot_L))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_R,t,ones(((nr_points*2)*nr_wb),1)*mean(theta_ddot_R))
%     subplot(3,1,2); plot(t,eta_ddot_R,t,ones(((nr_points*2)*nr_wb),1)*mean(eta_ddot_R))
%     subplot(3,1,3); plot(t,phi_ddot_R,t,ones(((nr_points*2)*nr_wb),1)*mean(phi_ddot_R))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL(1,:),t,ones(((nr_points*2)*nr_wb),1)*mean(wL(1,:)))
%     subplot(3,1,2); plot(t,wL(2,:),t,ones(((nr_points*2)*nr_wb),1)*mean(wL(2,:)))
%     subplot(3,1,3); plot(t,wL(3,:),t,ones(((nr_points*2)*nr_wb),1)*mean(wL(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wR(1,:),t,ones(((nr_points*2)*nr_wb),1)*mean(wR(1,:)))
%     subplot(3,1,2); plot(t,wR(2,:),t,ones(((nr_points*2)*nr_wb),1)*mean(wR(2,:)))
%     subplot(3,1,3); plot(t,wR(3,:),t,ones(((nr_points*2)*nr_wb),1)*mean(wR(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L(1,:),t,ones(((nr_points*2)*nr_wb),1)*mean(w_dot_L(1,:)))
%     subplot(3,1,2); plot(t,w_dot_L(2,:),t,ones(((nr_points*2)*nr_wb),1)*mean(w_dot_L(2,:)))
%     subplot(3,1,3); plot(t,w_dot_L(3,:),t,ones(((nr_points*2)*nr_wb),1)*mean(w_dot_L(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_R(1,:),t,ones(((nr_points*2)*nr_wb),1)*mean(w_dot_R(1,:)))
%     subplot(3,1,2); plot(t,w_dot_R(2,:),t,ones(((nr_points*2)*nr_wb),1)*mean(w_dot_R(2,:)))
%     subplot(3,1,3); plot(t,w_dot_R(3,:),t,ones(((nr_points*2)*nr_wb),1)*mean(w_dot_R(3,:)))
%     hold off
%     
    figure()
    hold on
    subplot(3,1,1); plot(t,wL_b(1,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(wL_b(1,:)))
    subplot(3,1,2); plot(t,wL_b(2,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(wL_b(2,:)))
    subplot(3,1,3); plot(t,wL_b(3,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(wL_b(3,:)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,wR_b(1,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(wR_b(1,:)))
    subplot(3,1,2); plot(t,wR_b(2,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(wR_b(2,:)))
    subplot(3,1,3); plot(t,wR_b(3,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(wR_b(3,:)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,w_dot_L_b(1,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(w_dot_L_b(1,:)))
    subplot(3,1,2); plot(t,w_dot_L_b(2,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(w_dot_L_b(2,:)))
    subplot(3,1,3); plot(t,w_dot_L_b(3,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(w_dot_L_b(3,:)))
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,w_dot_R_b(1,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(w_dot_R_b(1,:)))
    subplot(3,1,2); plot(t,w_dot_R_b(2,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(w_dot_R_b(2,:)))
    subplot(3,1,3); plot(t,w_dot_R_b(3,:),t,ones(((nr_points*2)*nr_wb+1),1)*mean(w_dot_R_b(3,:)))
    hold off
    
    wtL = [0; -1; 0];
    wtR = [0; 1; 0];
    
    xyz_L = zeros(3,(nr_points*nr_wb*2+1));
    xyz_R = zeros(3,(nr_points*nr_wb*2+1));
    
    for k = 1:(nr_points*nr_wb*2+1)
        xyz_L(:,k) = RL(:,:,k)'*wtL;
        xyz_R(:,k) = RR(:,:,k)'*wtR;
    end    
    
    figure()
    hold on
    subplot(3,1,1); plot(1:(nr_points*nr_wb*2+1),xyz_L(1,:),'r',1:(nr_points*nr_wb*2+1),xyz_R(1,:),'g')
    subplot(3,1,2); plot(1:(nr_points*nr_wb*2+1),xyz_L(2,:),'r',1:(nr_points*nr_wb*2+1),xyz_R(2,:),'g')
    subplot(3,1,3); plot(1:(nr_points*nr_wb*2+1),xyz_L(3,:),'r',1:(nr_points*nr_wb*2+1),xyz_R(3,:),'g')
    hold off
    
    
    
%     
%     pause
% 
    
%     wtL1 = [0.05; -1; 0];
%     wtR1 = [0.05; 1; 0];
%     wtL2 = [-0.05; -1; 0];
%     wtR2 = [-0.05; 1; 0];
%     wtL3 = [0.025; -1; 0];
%     wtR3 = [0.025; 1; 0];
%     
%     figure()
%     hold on
%     for i = 1:10:((nr_points*2)*nr_wb)
%         xyz1_L = RL(:,:,i)'*wtL1;
%         xyz2_L = RL(:,:,i)'*wtL2;
%         xyz1_R = RR(:,:,i)'*wtR1;
%         xyz2_R = RR(:,:,i)'*wtR2;
%         plot3([xyz1_L(1) xyz2_L(1)],[xyz1_L(2) xyz2_L(2)],[xyz1_L(3) xyz2_L(3)],'r')
%         plot3([xyz1_R(1) xyz2_R(1)],[xyz1_R(2) xyz2_R(2)],[xyz1_R(3) xyz2_R(3)],'g')
%     end
%     axis equal
%     hold off
%     
%     pause

%     wing_kin.R_strk     = R_strk;
%     wing_kin.vb         = [0;0;0];
%     wing_kin.wb         = [0;0;0];
%     wing_kin.ab         = [0;0;0];
%     wing_kin.w_dot_b    = [0;0;0];
%     wing_kin.RL         = RL;
%     wing_kin.RR         = RR;
%     wing_kin.wL         = wL;
%     wing_kin.wR         = wR;
%     wing_kin.w_dot_L    = w_dot_L;
%     wing_kin.w_dot_R    = w_dot_R;
%     wing_kin.wL_b       = wL_b;
%     wing_kin.wR_b       = wR_b;
%     wing_kin.w_dot_L_b  = w_dot_L_b;
%     wing_kin.w_dot_R_b  = w_dot_R_b;
    
    wingbeat.R_strk     = R_strk;
    
    FI_vel      = nan(3,nr_points*nr_wb*2+1);
    MI_vel      = nan(3,nr_points*nr_wb*2+1);
    
    for i = 1:(nr_points*nr_wb*2+1)
        
        kine.R_strk     = R_strk;
        kine.vb         = [0;0;0];
        kine.wb         = [0;0;0];
        kine.ab         = [0;0;0];
        kine.w_dot_b    = [0;0;0];
        kine.RL         = RL(:,:,i);
        kine.RR         = RR(:,:,i);
        kine.wL         = wL(:,i);
        kine.wR         = wR(:,i);
        kine.w_dot_L    = w_dot_L(:,i);
        kine.w_dot_R    = w_dot_R(:,i);
        kine.wL_b       = wL_b(:,i);
        kine.wR_b       = wR_b(:,i);
        kine.w_dot_L_b  = w_dot_L_b(:,i);
        kine.w_dot_R_b  = w_dot_R_b(:,i);
    
        [ FM_b, FM_strkpln ] = Inertia_instantaneous( kine, body_model, wing_model );
        
        FI_vel(:,i) = FM_b.F_I_vel;
        MI_vel(:,i) = FM_b.M_I_vel;
    
    end
    
    figure()
    hold on
    subplot(4,1,1); plot(t,FI_vel(1,:),t,ones(nr_points*nr_wb*2+1,1)*mean(FI_vel(1,1:(nr_points*2))))
    subplot(4,1,2); plot(t,FI_vel(2,:),t,ones(nr_points*nr_wb*2+1,1)*mean(FI_vel(2,1:(nr_points*2))))
    subplot(4,1,3); plot(t,FI_vel(3,:),t,ones(nr_points*nr_wb*2+1,1)*mean(FI_vel(3,1:(nr_points*2))))
    subplot(4,1,4); plot(t,(FI_vel(1,:).^2+FI_vel(2,:).^2+FI_vel(3,:).^2).^0.5)
    hold off
    
    figure()
    hold on
    subplot(3,1,1); plot(t,MI_vel(1,:),t,ones(nr_points*nr_wb*2+1,1)*mean(MI_vel(1,1:(nr_points*2))))
    subplot(3,1,2); plot(t,MI_vel(2,:),t,ones(nr_points*nr_wb*2+1,1)*mean(MI_vel(2,1:(nr_points*2))))
    subplot(3,1,3); plot(t,MI_vel(3,:),t,ones(nr_points*nr_wb*2+1,1)*mean(MI_vel(3,1:(nr_points*2))))
    hold off
    
    figure()
    hold on
    plot(MI_vel(2,:),'r')
    plot(ones(100,1)*mean(MI_vel(2,1:(2*nr_points))),'g')
    hold off
    
    mean(MI_vel(2,1:(nr_points+1)))/mean(FI_vel(3,1:(nr_points+1)))
    
%     % Initial conditions:
%     
%     state.dt            = dt;
%     state.id            = 1;
%     state.xyz           = zeros(3,1);
%     state.qb            = quat2mat(R_strk')';
%     state.Rb            = R_strk';
%     state.vb            = [0; 0; 0];
%     state.wb            = [0; 0; 0];
%     state.ab            = zeros(3,1);
%     state.w_dot_b       = zeros(3,1);
%     state.id_last       = nr_points*nr_wb;
%     
%     xyz                 = nan(3,nr_points*nr_wb);
%     qb                  = nan(4,nr_points*nr_wb);
%     Rb                  = nan(3,3,nr_points*nr_wb);
%     vb                  = nan(3,nr_points*nr_wb);
%     wb                  = nan(3,nr_points*nr_wb);
%     ab                  = nan(3,nr_points*nr_wb);
%     w_dot_b             = nan(3,nr_points*nr_wb);
%     
%     FI_vel              = nan(3,nr_points*nr_wb);
%     FI_acc              = nan(3,nr_points*nr_wb);
%     MI_vel              = nan(3,nr_points*nr_wb);
%     MI_acc              = nan(3,nr_points*nr_wb);
%     
%     for i = 1:(nr_wb*nr_points)
%         
%         i
%         
%         [ state_new ] = RK_updater_inertia_only( state, wing_kin, body_model, wing_model, wingbeat );
%         
%         state.id        = state_new.id;
%         state.xyz       = state_new.xyz;
%         state.qb        = state_new.qb;
%         state.Rb        = state_new.Rb;
%         state.vb        = state_new.vb;
%         state.wb        = state_new.wb;
%         state.ab        = state_new.ab;
%         state.w_dot_b   = state_new.w_dot_b;
%         
%         FI_vel(:,i)     = state_new.FI_vel;
%         MI_vel(:,i)     = state_new.MI_vel;
%         
%         xyz(:,i)        = state_new.xyz;
%         vb(:,i)         = state_new.vb;
%         wb(:,i)         = state_new.wb;
%         ab(:,i)         = state_new.ab;
%         w_dot_b(:,i)    = state_new.w_dot_b;
%     
%     end
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t(1:2:end),xyz(1,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(xyz(1,:)))
%     subplot(3,1,2); plot(t(1:2:end),xyz(2,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(xyz(2,:)))
%     subplot(3,1,3); plot(t(1:2:end),xyz(3,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(xyz(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t(1:2:end),vb(1,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(vb(1,:)))
%     subplot(3,1,2); plot(t(1:2:end),vb(2,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(vb(2,:)))
%     subplot(3,1,3); plot(t(1:2:end),vb(3,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(vb(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t(1:2:end),wb(1,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(wb(1,:)))
%     subplot(3,1,2); plot(t(1:2:end),wb(2,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(wb(2,:)))
%     subplot(3,1,3); plot(t(1:2:end),wb(3,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(wb(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t(1:2:end),FI_vel(1,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(FI_vel(1,:)))
%     subplot(3,1,2); plot(t(1:2:end),FI_vel(2,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(FI_vel(2,:)))
%     subplot(3,1,3); plot(t(1:2:end),FI_vel(3,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(FI_vel(3,:)))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t(1:2:end),MI_vel(1,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(MI_vel(1,:)))
%     subplot(3,1,2); plot(t(1:2:end),MI_vel(2,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(MI_vel(2,:)))
%     subplot(3,1,3); plot(t(1:2:end),MI_vel(3,:),t(1:2:end),ones(nr_points*nr_wb,1).*nanmean(MI_vel(3,:)))
%     hold off

end
