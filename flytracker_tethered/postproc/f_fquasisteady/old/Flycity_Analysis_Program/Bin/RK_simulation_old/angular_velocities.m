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
       
    d_theta_L_dt    = nan(1,nr_points);
    d_eta_L_dt      = nan(1,nr_points);
    d_phi_L_dt      = nan(1,nr_points);
    
    d_theta_R_dt    = nan(1,nr_points);
    d_eta_R_dt      = nan(1,nr_points);
    d_phi_R_dt      = nan(1,nr_points);
    
    d2_theta_L_dt2  = nan(1,nr_points);
    d2_eta_L_dt2    = nan(1,nr_points);
    d2_phi_L_dt2    = nan(1,nr_points);
    
    d2_theta_R_dt2  = nan(1,nr_points);
    d2_eta_R_dt2    = nan(1,nr_points);
    d2_phi_R_dt2    = nan(1,nr_points);
    
    for i = 1:nr_points
        
        if i == 1
        
            d_theta_L_dt(i)    = (theta_L(i+1)-theta_L(i))/dt;
            d_eta_L_dt(i)      = (eta_L(i+1)-eta_L(i))/dt;
            d_phi_L_dt(i)      = (phi_L(i+1)-phi_L(i))/dt;

            d_theta_R_dt(i)    = (theta_R(i+1)-theta_R(i))/dt;
            d_eta_R_dt(i)      = (eta_R(i+1)-eta_R(i))/dt;
            d_phi_R_dt(i)      = (phi_R(i+1)-phi_R(i))/dt;

            d2_theta_L_dt2(i)  = (theta_L(i+2)-2*theta_L(i+1)+theta_L(i))/(dt^2);
            d2_eta_L_dt2(i)    = (eta_L(i+2)-2*eta_L(i+1)+eta_L(i))/(dt^2);
            d2_phi_L_dt2(i)    = (phi_L(i+2)-2*phi_L(i+1)+phi_L(i))/(dt^2);

            d2_theta_R_dt2(i)  = (theta_R(i+2)-2*theta_R(i+1)+theta_R(i))/(dt^2);
            d2_eta_R_dt2(i)    = (eta_R(i+2)-2*eta_R(i+1)+eta_R(i))/(dt^2);
            d2_phi_R_dt2(i)    = (phi_R(i+2)-2*phi_R(i+1)+phi_R(i))/(dt^2);
        
        elseif i == nr_points

            d_theta_L_dt(i)    = (theta_L(i)-theta_L(i-1))/dt;
            d_eta_L_dt(i)      = (eta_L(i)-eta_L(i-1))/dt;
            d_phi_L_dt(i)      = (phi_L(i)-phi_L(i-1))/dt;

            d_theta_R_dt(i)    = (theta_R(i)-theta_R(i-1))/dt;
            d_eta_R_dt(i)      = (eta_R(i)-eta_R(i-1))/dt;
            d_phi_R_dt(i)      = (phi_R(i)-phi_R(i-1))/dt;

            d2_theta_L_dt2(i)  = (theta_L(i)-2*theta_L(i-1)+theta_L(i-2))/(dt^2);
            d2_eta_L_dt2(i)    = (eta_L(i)-2*eta_L(i-1)+eta_L(i-2))/(dt^2);
            d2_phi_L_dt2(i)    = (phi_L(i)-2*phi_L(i-1)+phi_L(i-2))/(dt^2);

            d2_theta_R_dt2(i)  = (theta_R(i)-2*theta_R(i-1)+theta_R(i-2))/(dt^2);
            d2_eta_R_dt2(i)    = (eta_R(i)-2*eta_R(i-1)+eta_R(i-2))/(dt^2);
            d2_phi_R_dt2(i)    = (phi_R(i)-2*phi_R(i-1)+phi_R(i-2))/(dt^2);
            
        else
            
            d_theta_L_dt(i)    = (theta_L(i+1)-theta_L(i-1))/(2*dt);
            d_eta_L_dt(i)      = (eta_L(i+1)-eta_L(i-1))/(2*dt);
            d_phi_L_dt(i)      = (phi_L(i+1)-phi_L(i-1))/(2*dt);

            d_theta_R_dt(i)    = (theta_R(i+1)-theta_R(i-1))/(2*dt);
            d_eta_R_dt(i)      = (eta_R(i+1)-eta_R(i-1))/(2*dt);
            d_phi_R_dt(i)      = (phi_R(i+1)-phi_R(i-1))/(2*dt);

            d2_theta_L_dt2(i)  = (theta_L(i+1)-2*theta_L(i)+theta_L(i-1))/(dt^2);
            d2_eta_L_dt2(i)    = (eta_L(i+1)-2*eta_L(i)+eta_L(i-1))/(dt^2);
            d2_phi_L_dt2(i)    = (phi_L(i+1)-2*phi_L(i)+phi_L(i-1))/(dt^2);

            d2_theta_R_dt2(i)  = (theta_R(i+1)-2*theta_R(i)+theta_R(i-1))/(dt^2);
            d2_eta_R_dt2(i)    = (eta_R(i+1)-2*eta_R(i)+eta_R(i-1))/(dt^2);
            d2_phi_R_dt2(i)    = (phi_R(i+1)-2*phi_R(i)+phi_R(i-1))/(dt^2);
            
        end
        
    end
    
    d_theta_L_dt    = d_theta_L_dt;
    d_eta_L_dt      = -d_eta_L_dt;
    d_phi_L_dt      = -d_phi_L_dt;
    
    d_theta_R_dt    = -d_theta_R_dt;
    d_eta_R_dt      = -d_eta_R_dt;
    d_phi_R_dt      = d_phi_R_dt;
    
    d2_theta_L_dt2  = d2_theta_L_dt2;
    d2_eta_L_dt2    = -d2_eta_L_dt2;
    d2_phi_L_dt2    = -d2_phi_L_dt2;
    
    d2_theta_R_dt2  = -d2_theta_R_dt2;
    d2_eta_R_dt2    = -d2_eta_R_dt2;
    d2_phi_R_dt2    = d2_phi_R_dt2;
    
    
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
        
        wL(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + R_180*[0; eta_dot_L(i); 0];
        wR(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + R_180*[0; eta_dot_R(i); 0];
        
        w_dot_L(:,i) = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + R_180*[0; eta_ddot_L(i); 0];
        w_dot_R(:,i) = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + R_180*[0; eta_ddot_R(i); 0];

        
        wL_b(:,i) = RL(:,:,i)'*wL(:,i);
        wR_b(:,i) = RR(:,:,i)'*wR(:,i);

        w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
        w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
        
        wL_strk(:,i) = R_strk*wL_b(:,i);
        wR_strk(:,i) = R_strk*wR_b(:,i);
        
        w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
        w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);
        
    end
    
    wL_b_old        = wL_b;
    wR_b_old        = wR_b;
    w_dot_L_b_old   = w_dot_L_b;
    w_dot_R_b_old   = w_dot_R_b;
    
    % Set the mean of wL_by and wR_by to zero
    
    wLbx_mean = mean(wL_b(1,1:(nr_points-1)));
    wRbx_mean = mean(wR_b(1,1:(nr_points-1)));
    wLby_mean = mean(wL_b(2,1:(nr_points-1)));
    wRby_mean = mean(wR_b(2,1:(nr_points-1)));
    wLbz_mean = mean(wL_b(3,1:(nr_points-1)));
    wRbz_mean = mean(wR_b(3,1:(nr_points-1)));
    w_dot_Lbx_mean = mean(w_dot_L_b(1,1:(nr_points-1)));
    w_dot_Rbx_mean = mean(w_dot_R_b(1,1:(nr_points-1)));
    w_dot_Lby_mean = mean(w_dot_L_b(2,1:(nr_points-1)));
    w_dot_Rby_mean = mean(w_dot_R_b(2,1:(nr_points-1)));
    w_dot_Lbz_mean = mean(w_dot_L_b(3,1:(nr_points-1)));
    w_dot_Rbz_mean = mean(w_dot_R_b(3,1:(nr_points-1)));
    
    for i = 1:nr_points
        
        wL_0 = RL(:,:,i)*[0; wLby_mean; 0];
        wR_0 = RR(:,:,i)*[0; wRby_mean; 0];
        w_dot_L_0 = RL(:,:,i)*[0; w_dot_Lby_mean; 0];
        w_dot_R_0 = RR(:,:,i)*[0; w_dot_Rby_mean; 0];
        
        wL(:,i) = wL(:,i)-wL_0;
        wR(:,i) = wR(:,i)-wR_0;
        w_dot_L(:,i) = w_dot_L(:,i)-w_dot_L_0;
        w_dot_R(:,i) = w_dot_R(:,i)-w_dot_R_0;
        
        wL_b(:,i) = RL(:,:,i)'*wL(:,i);
        wR_b(:,i) = RR(:,:,i)'*wR(:,i);
        w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
        w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
        
        wL_strk(:,i) = R_strk*wL_b(:,i);
        wR_strk(:,i) = R_strk*wR_b(:,i);
        w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
        w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);

%         wL_0 = RL(:,:,i)*[wLbx_mean; wLby_mean; wLbz_mean];
%         wR_0 = RR(:,:,i)*[wRbx_mean; wRby_mean; wRbz_mean];
%         w_dot_L_0 = RL(:,:,i)*[w_dot_Lbx_mean; w_dot_Lby_mean; w_dot_Lbz_mean];
%         w_dot_R_0 = RR(:,:,i)*[w_dot_Rbx_mean; w_dot_Rby_mean; w_dot_Rbz_mean];
%         
%         wL(:,i) = wL(:,i)-wL_0;
%         wR(:,i) = wR(:,i)-wR_0;
%         w_dot_L(:,i) = w_dot_L(:,i)-w_dot_L_0;
%         w_dot_R(:,i) = w_dot_R(:,i)-w_dot_R_0;
%         
%         wL_b(:,i) = RL(:,:,i)'*wL(:,i);
%         wR_b(:,i) = RR(:,:,i)'*wR(:,i);
%         w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
%         w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
%         
%         wL_strk(:,i) = R_strk*wL_b(:,i);
%         wR_strk(:,i) = R_strk*wR_b(:,i);
%         w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
%         w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);
        
    end
    
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
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_L,'r',t(m1),theta_L(m1),'o',t,ones(nr_points,1)*mean(theta_L),'k')
%     subplot(3,1,2); plot(t,eta_L,'r',t(m1),eta_L(m1),'o',t,ones(nr_points,1)*mean(eta_L),'k')
%     subplot(3,1,3); plot(t,phi_L,'r',t(m1),phi_L(m1),'o',t,ones(nr_points,1)*mean(phi_L),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_R,'r',t(m1),theta_R(m1),'o',t,ones(nr_points,1)*mean(theta_R),'k')
%     subplot(3,1,2); plot(t,eta_R,'r',t(m1),eta_R(m1),'o',t,ones(nr_points,1)*mean(eta_R),'k')
%     subplot(3,1,3); plot(t,phi_R,'r',t(m1),phi_R(m1),'o',t,ones(nr_points,1)*mean(phi_R),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_L,'b',t,d_theta_L_dt,'r',t(m1),theta_dot_L(m1),'o',t,ones(nr_points,1).*mean(theta_dot_L),'k')
%     subplot(3,1,2); plot(t,eta_dot_L,'b',t,d_eta_L_dt,'r',t(m1),eta_dot_L(m1),'o',t,ones(nr_points,1).*mean(eta_dot_L),'k')
%     subplot(3,1,3); plot(t,phi_dot_L,'b',t,d_phi_L_dt,'r',t(m1),phi_dot_L(m1),'o',t,ones(nr_points,1).*mean(phi_dot_L),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_R,'b',t,d_theta_R_dt,'r',t(m1),theta_dot_R(m1),'o',t,ones(nr_points,1).*mean(theta_dot_R),'k')
%     subplot(3,1,2); plot(t,eta_dot_R,'b',t,d_eta_R_dt,'r',t(m1),eta_dot_R(m1),'o',t,ones(nr_points,1).*mean(eta_dot_R),'k')
%     subplot(3,1,3); plot(t,phi_dot_R,'b',t,d_phi_R_dt,'r',t(m1),phi_dot_R(m1),'o',t,ones(nr_points,1).*mean(phi_dot_R),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_L,'b',t,d2_theta_L_dt2,'r',t(m1),theta_ddot_L(m1),'o',t,ones(nr_points,1).*mean(theta_ddot_L),'k')
%     subplot(3,1,2); plot(t,eta_ddot_L,'b',t,d2_eta_L_dt2,'r',t(m1),eta_ddot_L(m1),'o',t,ones(nr_points,1).*mean(eta_ddot_L),'k')
%     subplot(3,1,3); plot(t,phi_ddot_L,'b',t,d2_phi_L_dt2,'r',t(m1),phi_ddot_L(m1),'o',t,ones(nr_points,1).*mean(phi_ddot_L),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_R,'b',t,d2_theta_R_dt2,'r',t(m1),theta_ddot_R(m1),'o',t,ones(nr_points,1).*mean(theta_ddot_R),'k')
%     subplot(3,1,2); plot(t,eta_ddot_R,'b',t,d2_eta_R_dt2,'r',t(m1),eta_ddot_R(m1),'o',t,ones(nr_points,1).*mean(eta_ddot_R),'k')
%     subplot(3,1,3); plot(t,phi_ddot_R,'b',t,d2_phi_R_dt2,'r',t(m1),phi_ddot_R(m1),'o',t,ones(nr_points,1).*mean(phi_ddot_R),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL(1,:),'r',t,ones(nr_points,1).*mean(wL(1,:)),'k',t,wR(1,:),'b',t,ones(nr_points,1).*mean(wR(1,:)),'k')
%     subplot(3,1,2); plot(t,wL(2,:),'r',t,ones(nr_points,1).*mean(wL(2,:)),'k',t,wR(2,:),'b',t,ones(nr_points,1).*mean(wR(2,:)),'k')
%     subplot(3,1,3); plot(t,wL(3,:),'r',t,ones(nr_points,1).*mean(wL(3,:)),'k',t,wR(3,:),'b',t,ones(nr_points,1).*mean(wR(3,:)),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L(1,:),'r',t,ones(nr_points,1).*mean(w_dot_L(1,:)),'k',t,w_dot_R(1,:),'b',t,ones(nr_points,1).*mean(w_dot_R(1,:)),'k')
%     subplot(3,1,2); plot(t,w_dot_L(2,:),'r',t,ones(nr_points,1).*mean(w_dot_L(2,:)),'k',t,w_dot_R(2,:),'b',t,ones(nr_points,1).*mean(w_dot_R(2,:)),'k')
%     subplot(3,1,3); plot(t,w_dot_L(3,:),'r',t,ones(nr_points,1).*mean(w_dot_L(3,:)),'k',t,w_dot_R(3,:),'b',t,ones(nr_points,1).*mean(w_dot_R(3,:)),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL_b(1,:),'r',t,ones(nr_points,1).*mean(wL_b(1,:)),'k',t,wR_b(1,:),'b',t,ones(nr_points,1).*mean(wR_b(1,:)),'k',t,wL_b_old(1,:),'g')
%     subplot(3,1,2); plot(t,wL_b(2,:),'r',t,ones(nr_points,1).*mean(wL_b(2,:)),'k',t,wR_b(2,:),'b',t,ones(nr_points,1).*mean(wR_b(2,:)),'k',t,wL_b_old(2,:),'g')
%     subplot(3,1,3); plot(t,wL_b(3,:),'r',t,ones(nr_points,1).*mean(wL_b(3,:)),'k',t,wR_b(3,:),'b',t,ones(nr_points,1).*mean(wR_b(3,:)),'k',t,wL_b_old(3,:),'g')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L_b(1,:),'r',t,ones(nr_points,1).*mean(w_dot_L_b(1,:)),'k',t,w_dot_R_b(1,:),'b',t,ones(nr_points,1).*mean(w_dot_R_b(1,:)),'k')
%     subplot(3,1,2); plot(t,w_dot_L_b(2,:),'r',t,ones(nr_points,1).*mean(w_dot_L_b(2,:)),'k',t,w_dot_R_b(2,:),'b',t,ones(nr_points,1).*mean(w_dot_R_b(2,:)),'k')
%     subplot(3,1,3); plot(t,w_dot_L_b(3,:),'r',t,ones(nr_points,1).*mean(w_dot_L_b(3,:)),'k',t,w_dot_R_b(3,:),'b',t,ones(nr_points,1).*mean(w_dot_R_b(3,:)),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL_strk(1,:),'r',t,ones(nr_points,1).*mean(wL_strk(1,:)),'k',t,wR_strk(1,:),'b',t,ones(nr_points,1).*mean(wR_strk(1,:)),'k')
%     subplot(3,1,2); plot(t,wL_strk(2,:),'r',t,ones(nr_points,1).*mean(wL_strk(2,:)),'k',t,wR_strk(2,:),'b',t,ones(nr_points,1).*mean(wR_strk(2,:)),'k')
%     subplot(3,1,3); plot(t,wL_strk(3,:),'r',t,ones(nr_points,1).*mean(wL_strk(3,:)),'k',t,wR_strk(3,:),'b',t,ones(nr_points,1).*mean(wR_strk(3,:)),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L_strk(1,:),'r',t,ones(nr_points,1).*mean(w_dot_L_strk(1,:)),'k',t,w_dot_R_strk(1,:),'b',t,ones(nr_points,1).*mean(w_dot_R_strk(1,:)),'k')
%     subplot(3,1,2); plot(t,w_dot_L_strk(2,:),'r',t,ones(nr_points,1).*mean(w_dot_L_strk(2,:)),'k',t,w_dot_R_strk(2,:),'b',t,ones(nr_points,1).*mean(w_dot_R_strk(2,:)),'k')
%     subplot(3,1,3); plot(t,w_dot_L_strk(3,:),'r',t,ones(nr_points,1).*mean(w_dot_L_strk(3,:)),'k',t,w_dot_R_strk(3,:),'b',t,ones(nr_points,1).*mean(w_dot_R_strk(3,:)),'k')
%     hold off
%     
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
%     for i = 1:round(nr_points/50):(nr_points)
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
%     

end

