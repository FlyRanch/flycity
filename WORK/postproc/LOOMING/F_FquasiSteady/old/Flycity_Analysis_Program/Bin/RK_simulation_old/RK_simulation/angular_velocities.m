function [ kine ] = angular_velocities( a_fit )


    % Compute the wing kinematics:
    
    a_theta_L   = a_fit.a_theta_L;
    a_eta_L     = a_fit.a_eta_L;
    a_phi_L     = a_fit.a_phi_L;
    a_theta_R   = a_fit.a_theta_R;
    a_eta_R     = a_fit.a_eta_R;
    a_phi_R     = a_fit.a_phi_R;
    f           = a_fit.f;
    down_up     = a_fit.down_up;
    nr_points   = a_fit.nr_points;
    R_strk      = a_fit.R_strk;

    n_pol_theta     = (length(a_theta_L)-2)/2;
    n_pol_eta       = (length(a_eta_L)-2)/2;
    n_pol_phi       = (length(a_phi_L)-2)/2;
    
    nr_points_down  = round(down_up*nr_points);
    
    dt = (1/f)/(nr_points-1);
    
    [ t, X_theta ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up, nr_points, 0, dt*(nr_points-1), 0 );
    [ ~, X_eta ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up, nr_points, 0, dt*(nr_points-1), 0 );
    [ ~, X_phi ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up, nr_points, 0, dt*(nr_points-1), 0 );
    
    [ ~, X_theta_dot ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up, nr_points, 0, dt*(nr_points-1), 1 );
    [ ~, X_eta_dot ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up, nr_points, 0, dt*(nr_points-1), 1 );
    [ ~, X_phi_dot ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up, nr_points, 0, dt*(nr_points-1), 1 );
    
    [ ~, X_theta_ddot ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up, nr_points, 0, dt*(nr_points-1), 2 );
    [ ~, X_eta_ddot ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up, nr_points, 0, dt*(nr_points-1), 2 );
    [ ~, X_phi_ddot ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up, nr_points, 0, dt*(nr_points-1), 2 );    

    
    m1 = nr_points_down;
    m2 = nr_points-nr_points_down+1;
    
    theta_L         = X_theta*a_theta_L;
    eta_L           = X_eta*a_eta_L;
    phi_L           = X_phi*a_phi_L;
    
    theta_R         = X_theta*a_theta_R;
    eta_R           = X_eta*a_eta_R;
    phi_R           = X_phi*a_phi_R;
    
    theta_dot_L_t   = X_theta_dot*a_theta_L;
    eta_dot_L_t     = X_eta_dot*a_eta_L;
    phi_dot_L_t     = X_phi_dot*a_phi_L;
    
    theta_dot_R_t   = X_theta_dot*a_theta_R;
    eta_dot_R_t     = X_eta_dot*a_eta_R;
    phi_dot_R_t     = X_phi_dot*a_phi_R;
    
    theta_ddot_L_t  = X_theta_ddot*a_theta_L;
    eta_ddot_L_t    = X_eta_ddot*a_eta_L;
    phi_ddot_L_t    = X_phi_ddot*a_phi_L;
    
    theta_ddot_R_t  = X_theta_ddot*a_theta_R;
    eta_ddot_R_t    = X_eta_ddot*a_eta_R;
    phi_ddot_R_t    = X_phi_ddot*a_phi_R;
    
%     theta_dot_L           = (theta_dot_L_t-mean(theta_dot_L_t))*(2*f);
%     eta_dot_L             = -(eta_dot_L_t-mean(eta_dot_L_t))*(2*f);
%     phi_dot_L             = -(phi_dot_L_t-mean(phi_dot_L_t))*(2*f);
%     
%     theta_dot_R           = -(theta_dot_R_t-mean(theta_dot_R_t))*(2*f);
%     eta_dot_R             = -(eta_dot_R_t-mean(eta_dot_R_t))*(2*f);
%     phi_dot_R             = (phi_dot_R_t-mean(phi_dot_R_t))*(2*f);
%     
%     theta_ddot_L           = (theta_ddot_L_t-mean(theta_ddot_L_t))*(2*f)^2;
%     eta_ddot_L             = -(eta_ddot_L_t-mean(eta_ddot_L_t))*(2*f)^2;
%     phi_ddot_L             = -(phi_ddot_L_t-mean(phi_ddot_L_t))*(2*f)^2;
%     
%     theta_ddot_R           = -(theta_ddot_R_t-mean(theta_ddot_R_t))*(2*f)^2;
%     eta_ddot_R             = -(eta_ddot_R_t-mean(eta_ddot_R_t))*(2*f)^2;
%     phi_ddot_R             = (phi_ddot_R_t-mean(phi_ddot_R_t))*(2*f)^2;


    theta_dot_L           = (theta_dot_L_t)*(2*f);
    eta_dot_L             = -(eta_dot_L_t)*(2*f);
    phi_dot_L             = -(phi_dot_L_t)*(2*f);
    
    theta_dot_R           = -(theta_dot_R_t)*(2*f);
    eta_dot_R             = -(eta_dot_R_t)*(2*f);
    phi_dot_R             = (phi_dot_R_t)*(2*f);
    
    theta_ddot_L           = (theta_ddot_L_t)*(2*f)^2;
    eta_ddot_L             = -(eta_ddot_L_t)*(2*f)^2;
    phi_ddot_L             = -(phi_ddot_L_t)*(2*f)^2;
    
    theta_ddot_R           = -(theta_ddot_R_t)*(2*f)^2;
    eta_ddot_R             = -(eta_ddot_R_t)*(2*f)^2;
    phi_ddot_R             = (phi_ddot_R_t)*(2*f)^2;

    
    R_phi_L     = nan(3,3,nr_points);
    R_theta_L   = nan(3,3,nr_points);
    R_eta_L     = nan(3,3,nr_points);
    
    R_phi_R     = nan(3,3,nr_points);
    R_theta_R   = nan(3,3,nr_points);
    R_eta_R     = nan(3,3,nr_points);
    
    RL = nan(3,3,nr_points);
    RR = nan(3,3,nr_points);
    
    wL = nan(3,nr_points);
    wR = nan(3,nr_points);
    
    w_dot_L = nan(3,nr_points);
    w_dot_R = nan(3,nr_points);
    
    wL_b = nan(3,nr_points);
    wR_b = nan(3,nr_points);
    
    w_dot_L_b = nan(3,nr_points);
    w_dot_R_b = nan(3,nr_points);
    
    wL_strk = nan(3,nr_points);
    wR_strk = nan(3,nr_points);
    
    w_dot_L_strk = nan(3,nr_points);
    w_dot_R_strk = nan(3,nr_points);
       
    for i = 1:nr_points
        
        % Matrix to switch orientation wing:
        
        R_180              = [ -1  0  0; ...
                                0  1  0; ...
                                0  0 -1];
                           
        % Wing kinematics rotation matrices:
        
        R_phi_L(:,:,i)     = [ cos(phi_L(i)) -sin(phi_L(i)) 0; ...
                               sin(phi_L(i))  cos(phi_L(i)) 0; ...
                               0              0             1];

        R_theta_L(:,:,i)   = [ 1  0               0              ; ...
                               0  cos(theta_L(i)) sin(theta_L(i)); ...
                               0 -sin(theta_L(i)) cos(theta_L(i))];

        R_eta_L(:,:,i)     = [  cos(eta_L(i)) 0 sin(eta_L(i)); ...
                                0             1 0            ; ...
                               -sin(eta_L(i)) 0 cos(eta_L(i))];
    
        R_phi_R(:,:,i)     = [ cos(phi_R(i)) sin(phi_R(i)) 0; ...
                              -sin(phi_R(i)) cos(phi_R(i)) 0; ...
                               0             0             1];

        R_theta_R(:,:,i)   = [ 1 0                0              ; ...
                               0 cos(theta_R(i)) -sin(theta_R(i)); ...
                               0 sin(theta_R(i))  cos(theta_R(i))];

        R_eta_R(:,:,i)     = [  cos(eta_R(i)) 0 sin(eta_R(i)); ...
                                0             1 0            ; ...
                               -sin(eta_R(i)) 0 cos(eta_R(i))];
                    
        RL(:,:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*R_phi_L(:,:,i)*R_strk;
        RR(:,:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*R_phi_R(:,:,i)*R_strk;
        
        RL(:,:,i) = RL(:,:,i)./norm(RL(:,:,i));
        RR(:,:,i) = RR(:,:,i)./norm(RR(:,:,i));

        wL(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + R_180*[0; eta_dot_L(i); 0];
        wR(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + R_180*[0; eta_dot_R(i); 0];
           
        w_dot_L(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + R_180*[0; eta_ddot_L(i); 0];
        w_dot_R(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + R_180*[0; eta_ddot_R(i); 0];

        
    end
    

    
%     wL_b(1,:) = wL_b(1,:)-mean(wL_b(1,1:(end-1)));
%     wL_b(2,:) = wL_b(2,:)-mean(wL_b(2,1:(end-1)));
%     wL_b(3,:) = wL_b(3,:)-mean(wL_b(3,1:(end-1)));
%     wR_b(1,:) = wR_b(1,:)-mean(wR_b(1,1:(end-1)));
%     wR_b(2,:) = wR_b(2,:)-mean(wR_b(2,1:(end-1)));
%     wR_b(3,:) = wR_b(3,:)-mean(wR_b(3,1:(end-1)));
%     w_dot_L_b(1,:) = w_dot_L_b(1,:)-mean(w_dot_L_b(1,1:(end-1)));
%     w_dot_L_b(2,:) = w_dot_L_b(2,:)-mean(w_dot_L_b(2,1:(end-1)));
%     w_dot_L_b(3,:) = w_dot_L_b(3,:)-mean(w_dot_L_b(3,1:(end-1)));
%     w_dot_R_b(1,:) = w_dot_R_b(1,:)-mean(w_dot_R_b(1,1:(end-1)));
%     w_dot_R_b(2,:) = w_dot_R_b(2,:)-mean(w_dot_R_b(2,1:(end-1)));
%     w_dot_R_b(3,:) = w_dot_R_b(3,:)-mean(w_dot_R_b(3,1:(end-1)));
    
    
    kine.dt             = dt;
    kine.t              = t;
    kine.theta_L        = theta_L;
    kine.eta_L          = eta_L;
    kine.phi_L          = phi_L;
    kine.theta_R        = theta_R;
    kine.eta_R          = eta_R;
    kine.phi_R          = phi_R;
    kine.theta_dot_L    = theta_dot_L;
    kine.eta_dot_L      = eta_dot_L;
    kine.phi_dot_L      = phi_dot_L;
    kine.theta_dot_R    = theta_dot_R;
    kine.eta_dot_R      = eta_dot_R;
    kine.phi_dot_R      = phi_dot_R;
    kine.theta_ddot_L   = theta_ddot_L;
    kine.eta_ddot_L     = eta_ddot_L;
    kine.phi_ddot_L     = phi_ddot_L;
    kine.theta_ddot_R   = theta_ddot_R;
    kine.eta_ddot_R     = eta_ddot_R;
    kine.phi_ddot_R     = phi_ddot_R;
    kine.RL             = RL;
    kine.RR             = RR;
    kine.wL             = wL;
    kine.wR             = wR;
    kine.w_dot_L        = w_dot_L;
    kine.w_dot_R        = w_dot_R;
    kine.wL_b           = wL_b;
    kine.wR_b           = wR_b;
    kine.w_dot_L_b      = w_dot_L_b;
    kine.w_dot_R_b      = w_dot_R_b;
    kine.wL_strk        = wL_strk;
    kine.wR_strk        = wR_strk;
    kine.w_dot_L_strk   = w_dot_L_strk;
    kine.w_dot_R_strk   = w_dot_R_strk;


end

