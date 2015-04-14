function [ kine ] = angular_velocities_polynomial( a_fit )


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
    
    theta_L_t       = X_theta*a_theta_L;
    eta_L_t         = -X_eta*a_eta_L;
    phi_L_t         = -X_phi*a_phi_L;
    
    theta_R_t       = -X_theta*a_theta_R;
    eta_R_t         = -X_eta*a_eta_R;
    phi_R_t         = X_phi*a_phi_R;
    
    theta_dot_L_t   = X_theta_dot*a_theta_L*(2*f);
    eta_dot_L_t     = -X_eta_dot*a_eta_L*(2*f);
    phi_dot_L_t     = -X_phi_dot*a_phi_L*(2*f);
    
    theta_dot_R_t   = -X_theta_dot*a_theta_R*(2*f);
    eta_dot_R_t     = -X_eta_dot*a_eta_R*(2*f);
    phi_dot_R_t     = X_phi_dot*a_phi_R*(2*f);
    
    theta_ddot_L_t  = X_theta_ddot*a_theta_L*(2*f)^2;
    eta_ddot_L_t    = -X_eta_ddot*a_eta_L*(2*f)^2;
    phi_ddot_L_t    = -X_phi_ddot*a_phi_L*(2*f)^2;
    
    theta_ddot_R_t  = -X_theta_ddot*a_theta_R*(2*f)^2;
    eta_ddot_R_t    = -X_eta_ddot*a_eta_R*(2*f)^2;
    phi_ddot_R_t    = X_phi_ddot*a_phi_R*(2*f)^2;



%     theta_L_t       = -X_theta*a_theta_L;
%     eta_L_t         = -X_eta*a_eta_L;
%     phi_L_t         = X_phi*a_phi_L;
%     
%     theta_R_t       = X_theta*a_theta_R;
%     eta_R_t         = -X_eta*a_eta_R;
%     phi_R_t         = -X_phi*a_phi_R;
%     
%     theta_dot_L_t   = -X_theta_dot*a_theta_L*(2*f);
%     eta_dot_L_t     = -X_eta_dot*a_eta_L*(2*f);
%     phi_dot_L_t     = X_phi_dot*a_phi_L*(2*f);
%     
%     theta_dot_R_t   = X_theta_dot*a_theta_R*(2*f);
%     eta_dot_R_t     = -X_eta_dot*a_eta_R*(2*f);
%     phi_dot_R_t     = -X_phi_dot*a_phi_R*(2*f);
%     
%     theta_ddot_L_t  = -X_theta_ddot*a_theta_L*(2*f)^2;
%     eta_ddot_L_t    = -X_eta_ddot*a_eta_L*(2*f)^2;
%     phi_ddot_L_t    = X_phi_ddot*a_phi_L*(2*f)^2;
%     
%     theta_ddot_R_t  = X_theta_ddot*a_theta_R*(2*f)^2;
%     eta_ddot_R_t    = -X_eta_ddot*a_eta_R*(2*f)^2;
%     phi_ddot_R_t    = -X_phi_ddot*a_phi_R*(2*f)^2;
    
%     d_theta_dt_L    = gradient(theta_L_t)./dt;
%     d_eta_dt_L      = gradient(eta_L_t)./dt;
%     d_phi_dt_L      = gradient(phi_L_t)./dt;
%     
%     d_theta_dt_R    = gradient(theta_R_t)./dt;
%     d_eta_dt_R      = gradient(eta_R_t)./dt;
%     d_phi_dt_R      = gradient(phi_R_t)./dt;
%     
%     d2_theta_dt2_L  = gradient(d_theta_dt_L)./dt;
%     d2_eta_dt2_L    = gradient(d_eta_dt_L)./dt;
%     d2_phi_dt2_L    = gradient(d_phi_dt_L)./dt;
%     
%     d2_theta_dt2_R  = gradient(d_theta_dt_R)./dt;
%     d2_eta_dt2_R    = gradient(d_eta_dt_R)./dt;
%     d2_phi_dt2_R    = gradient(d_phi_dt_R)./dt;
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_L_t,'r',t,theta_R_t,'b')
%     subplot(3,1,2); plot(t,eta_L_t,'r',t,eta_R_t,'b')
%     subplot(3,1,3); plot(t,phi_L_t,'r',t,phi_R_t,'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_L_t,'r',t,theta_dot_R_t,'b',t,d_theta_dt_L,'k',t,d_theta_dt_R,'k')
%     subplot(3,1,2); plot(t,eta_dot_L_t,'r',t,eta_dot_R_t,'b',t,d_eta_dt_L,'k',t,d_eta_dt_R,'k')
%     subplot(3,1,3); plot(t,phi_dot_L_t,'r',t,phi_dot_R_t,'b',t,d_phi_dt_L,'k',t,d_phi_dt_R,'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_L_t,'r',t,theta_ddot_R_t,'b',t,d2_theta_dt2_L,'k',t,d2_theta_dt2_R,'k')
%     subplot(3,1,2); plot(t,eta_ddot_L_t,'r',t,eta_ddot_R_t,'b',t,d2_eta_dt2_L,'k',t,d2_eta_dt2_R,'k')
%     subplot(3,1,3); plot(t,phi_ddot_L_t,'r',t,phi_ddot_R_t,'b',t,d2_phi_dt2_L,'k',t,d2_phi_dt2_R,'k')
%     hold off
%     
%     pause
    
%     b1_theta_L      = mean(theta_ddot_L_t);
%     b1_eta_L        = mean(eta_ddot_L_t);
%     b1_phi_L        = mean(phi_ddot_L_t);
%     
%     b1_theta_R      = mean(theta_ddot_R_t);
%     b1_eta_R        = mean(eta_ddot_R_t);
%     b1_phi_R        = mean(phi_ddot_R_t);
%     
%     theta_ddot_L    = theta_ddot_L_t-b1_theta_L;
%     eta_ddot_L      = eta_ddot_L_t-b1_eta_L;
%     phi_ddot_L      = phi_ddot_L_t-b1_phi_L;
%     
%     theta_ddot_R    = theta_ddot_R_t-b1_theta_R;
%     eta_ddot_R      = eta_ddot_R_t-b1_eta_R;
%     phi_ddot_R      = phi_ddot_R_t-b1_phi_R;
%     
%     theta_dot_L_t   = theta_dot_L_t-(1/(2*f))*b1_theta_L;
%     eta_dot_L_t     = eta_dot_L_t-(1/(2*f))*b1_eta_L;
%     phi_dot_L_t     = phi_dot_L_t-(1/(2*f))*b1_phi_L;
%     
%     theta_dot_R_t   = theta_dot_R_t-(1/(2*f))*b1_theta_R;
%     eta_dot_R_t     = eta_dot_R_t-(1/(2*f))*b1_eta_R;
%     phi_dot_R_t     = phi_dot_R_t-(1/(2*f))*b1_phi_R;
%     
%     b2_theta_L      = mean(theta_dot_L_t);
%     b2_eta_L        = mean(eta_dot_L_t);
%     b2_phi_L        = mean(phi_dot_L_t);
%     
%     b2_theta_R      = mean(theta_dot_R_t);
%     b2_eta_R        = mean(eta_dot_R_t);
%     b2_phi_R        = mean(phi_dot_R_t);
%     
%     theta_dot_L     = theta_dot_L_t-b2_theta_L;
%     eta_dot_L       = eta_dot_L_t-b2_eta_L;
%     phi_dot_L       = phi_dot_L_t-b2_phi_L;
%     
%     theta_dot_R     = theta_dot_R_t-b2_theta_R;
%     eta_dot_R       = eta_dot_R_t-b2_eta_R;
%     phi_dot_R       = phi_dot_R_t-b2_phi_R;
%     
%     theta_L         = theta_L_t-(1/(2*f))*b2_theta_L;
%     eta_L           = eta_L_t-(1/(2*f))*b2_eta_L;
%     phi_L           = phi_L_t-(1/(2*f))*b2_phi_L;
%     
%     theta_R         = theta_R_t-(1/(2*f))*b2_theta_R;
%     eta_R           = eta_R_t-(1/(2*f))*b2_eta_R;
%     phi_R           = phi_R_t-(1/(2*f))*b2_phi_R;

    theta_L         = theta_L_t;
    eta_L           = eta_L_t;
    phi_L           = phi_L_t;
    
    theta_R         = theta_R_t;
    eta_R           = eta_R_t;
    phi_R           = phi_R_t;
    
    theta_dot_L     = theta_dot_L_t;
    eta_dot_L       = eta_dot_L_t;
    phi_dot_L       = phi_dot_L_t;
    
    theta_dot_R     = theta_dot_R_t;
    eta_dot_R       = eta_dot_R_t;
    phi_dot_R       = phi_dot_R_t;
    
    theta_ddot_L    = theta_ddot_L_t;
    eta_ddot_L      = eta_ddot_L_t;
    phi_ddot_L      = phi_ddot_L_t;
    
    theta_ddot_R    = theta_ddot_R_t;
    eta_ddot_R      = eta_ddot_R_t;
    phi_ddot_R      = phi_ddot_R_t;
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_L,'r',t,theta_R,'b',t,ones(nr_points,1)*mean(theta_L),'k')
%     subplot(3,1,2); plot(t,eta_L,'r',t,eta_R,'b',t,ones(nr_points,1)*mean(eta_L),'k')
%     subplot(3,1,3); plot(t,phi_L,'r',t,phi_R,'b',t,ones(nr_points,1)*mean(phi_L),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_L,'r',t,theta_dot_R,'b',t,ones(nr_points,1)*mean(theta_dot_L),'k')
%     subplot(3,1,2); plot(t,eta_dot_L,'r',t,eta_dot_R,'b',t,ones(nr_points,1)*mean(eta_dot_L),'k')
%     subplot(3,1,3); plot(t,phi_dot_L,'r',t,phi_dot_R,'b',t,ones(nr_points,1)*mean(phi_dot_L),'k')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_L,'r',t,theta_ddot_R,'b',t,ones(nr_points,1)*mean(theta_ddot_L),'k')
%     subplot(3,1,2); plot(t,eta_ddot_L,'r',t,eta_ddot_R,'b',t,ones(nr_points,1)*mean(eta_ddot_L),'k')
%     subplot(3,1,3); plot(t,phi_ddot_L,'r',t,phi_ddot_R,'b',t,ones(nr_points,1)*mean(phi_ddot_L),'k')
%     hold off
%     
% %     pause
    
    R_phi_L         = nan(3,3,nr_points);
    R_theta_L       = nan(3,3,nr_points);
    R_eta_L         = nan(3,3,nr_points);
    
    R_phi_R         = nan(3,3,nr_points);
    R_theta_R       = nan(3,3,nr_points);
    R_eta_R         = nan(3,3,nr_points);
    
    R_dot_theta_L   = nan(3,3,nr_points);
    R_dot_eta_L     = nan(3,3,nr_points);
    
    R_dot_theta_R   = nan(3,3,nr_points);
    R_dot_eta_R     = nan(3,3,nr_points);
    
    RL = nan(3,3,nr_points);
    RR = nan(3,3,nr_points);
    
    qL = nan(4,nr_points);
    qR = nan(4,nr_points);
    
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
        
        R_phi_L(:,:,i)     = [ cos(phi_L(i))  sin(phi_L(i)) 0; ...
                              -sin(phi_L(i))  cos(phi_L(i)) 0; ...
                               0              0             1];

        R_theta_L(:,:,i)   = [ 1  0               0              ; ...
                               0  cos(theta_L(i)) sin(theta_L(i)); ...
                               0 -sin(theta_L(i)) cos(theta_L(i))];

        R_eta_L(:,:,i)     = [  cos(eta_L(i)) 0 -sin(eta_L(i)); ...
                                0             1  0            ; ...
                                sin(eta_L(i)) 0  cos(eta_L(i))];
    
        R_phi_R(:,:,i)     = [ cos(phi_R(i))  sin(phi_R(i)) 0; ...
                              -sin(phi_R(i))  cos(phi_R(i)) 0; ...
                               0              0             1];

        R_theta_R(:,:,i)   = [ 1  0               0              ; ...
                               0  cos(theta_R(i)) sin(theta_R(i)); ...
                               0 -sin(theta_R(i)) cos(theta_R(i))];

        R_eta_R(:,:,i)     = [  cos(eta_R(i)) 0 -sin(eta_R(i)); ...
                                0             1  0            ; ...
                                sin(eta_R(i)) 0  cos(eta_R(i))];
                           
        R_dot_eta_L(:,:,i) = [ -eta_dot_L(i)*sin(eta_L(i)) 0 -eta_dot_L(i)*cos(eta_L(i)); ...
                                0                          0  0                         ; ...
                                eta_dot_L(i)*cos(eta_L(i)) 0 -eta_dot_L(i)*sin(eta_L(i))];
                           
        R_dot_theta_L(:,:,i) = [ 0  0                               0                             ; ...
                                 0 -theta_dot_L(i)*sin(theta_L(i))  theta_dot_L(i)*cos(theta_L(i)); ...
                                 0 -theta_dot_L(i)*cos(theta_L(i)) -theta_dot_L(i)*sin(theta_L(i))];
                           
        R_dot_eta_R(:,:,i) = [ -eta_dot_R(i)*sin(eta_R(i)) 0 -eta_dot_R(i)*cos(eta_R(i)); ...
                                0                          0  0                         ; ...
                                eta_dot_R(i)*cos(eta_R(i)) 0 -eta_dot_R(i)*sin(eta_R(i))];
                           
        R_dot_theta_R(:,:,i) = [ 0  0                               0                             ; ...
                                 0 -theta_dot_R(i)*sin(theta_R(i))  theta_dot_R(i)*cos(theta_R(i)); ...
                                 0 -theta_dot_R(i)*cos(theta_R(i)) -theta_dot_R(i)*sin(theta_R(i))];

        RL(:,:,i)    = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*R_phi_L(:,:,i)*R_strk;
        RR(:,:,i)    = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*R_phi_R(:,:,i)*R_strk;
        
        RL(:,:,i)    = RL(:,:,i)./norm(RL(:,:,i));
        RR(:,:,i)    = RR(:,:,i)./norm(RR(:,:,i));
        
        qL_t         = quat2mat(RL(:,:,i));
        qR_t         = quat2mat(RR(:,:,i));
        
        qL(:,i)      = qL_t/norm(qL_t);
        qR(:,i)      = qR_t/norm(qR_t);

        wL(:,i)      = R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + R_180*R_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + R_180*[0; eta_dot_L(i); 0];
        wR(:,i)      = R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + R_180*R_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + R_180*[0; eta_dot_R(i); 0];
        
        w_dot_L(:,i) = R_180*R_dot_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + ...
                       R_180*R_eta_L(:,:,i)*R_dot_theta_L(:,:,i)*[0; 0; phi_dot_L(i)] + ...
                       R_180*R_eta_L(:,:,i)*R_theta_L(:,:,i)*[0; 0; phi_ddot_L(i)] + ...
                       R_180*R_dot_eta_L(:,:,i)*[theta_dot_L(i); 0; 0] + ...
                       R_180*R_eta_L(:,:,i)*[theta_ddot_L(i); 0; 0] + ...
                       R_180*[0; eta_ddot_L(i); 0];
        
        w_dot_R(:,i) = R_180*R_dot_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + ...
                       R_180*R_eta_R(:,:,i)*R_dot_theta_R(:,:,i)*[0; 0; phi_dot_R(i)] + ...
                       R_180*R_eta_R(:,:,i)*R_theta_R(:,:,i)*[0; 0; phi_ddot_R(i)] + ...
                       R_180*R_dot_eta_R(:,:,i)*[theta_dot_R(i); 0; 0] + ...
                       R_180*R_eta_R(:,:,i)*[theta_ddot_R(i); 0; 0] + ...
                       R_180*[0; eta_ddot_R(i); 0];
                   
        wL_b(:,i) = RL(:,:,i)'*wL(:,i);
        wR_b(:,i) = RR(:,:,i)'*wR(:,i);
        
        w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
        w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
        
        wL_strk(:,i) = R_strk*wL_b(:,i);
        wR_strk(:,i) = R_strk*wR_b(:,i);
        
        w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
        w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);


        
    end
    
%     wL_mean             = [mean(wL(1,1:(nr_points-1))); mean(wL(2,1:(nr_points-1))); mean(wL(3,1:(nr_points-1)))];
%     wR_mean             = [mean(wR(1,1:(nr_points-1))); mean(wR(2,1:(nr_points-1))); mean(wR(3,1:(nr_points-1)))];
%     w_dot_L_mean        = [mean(w_dot_L(1,1:(nr_points-1))); mean(w_dot_L(2,1:(nr_points-1))); mean(w_dot_L(3,1:(nr_points-1)))];
%     w_dot_R_mean        = [mean(w_dot_R(1,1:(nr_points-1))); mean(w_dot_R(2,1:(nr_points-1))); mean(w_dot_R(3,1:(nr_points-1)))];
%         
%     for i = 1:nr_points
% 
%         w_dot_L(:,i)        = w_dot_L(:,i)-w_dot_L_mean;
%         w_dot_R(:,i)        = w_dot_R(:,i)-w_dot_R_mean;
%         
%         wL(2,i)             = wL(2,i)-wL_mean(2);
%         wR(2,i)             = wR(2,i)-wR_mean(2);
%         
%     end
%     
%     for i = 1:nr_points
%         
%         wL_b(:,i) = RL(:,:,i)'*wL(:,i);
%         wR_b(:,i) = RR(:,:,i)'*wR(:,i);
%         
%         w_dot_L_b(:,i) = RL(:,:,i)'*w_dot_L(:,i);
%         w_dot_R_b(:,i) = RR(:,:,i)'*w_dot_R(:,i);
%         
%         wL_strk(:,i) = R_strk*wL_b(:,i);
%         wR_strk(:,i) = R_strk*wR_b(:,i);
%         
%         w_dot_L_strk(:,i) = R_strk*w_dot_L_b(:,i);
%         w_dot_R_strk(:,i) = R_strk*w_dot_R_b(:,i);        
%         
%     end
%     
%     w_dot_L_b_mean      = [mean(w_dot_L_b(1,1:(nr_points-1))); mean(w_dot_L_b(2,1:(nr_points-1))); mean(w_dot_L_b(3,1:(nr_points-1)))];
%     w_dot_R_b_mean      = [mean(w_dot_R_b(1,1:(nr_points-1))); mean(w_dot_R_b(2,1:(nr_points-1))); mean(w_dot_R_b(3,1:(nr_points-1)))];
%     w_dot_L_strk_mean   = [mean(w_dot_L_strk(1,1:(nr_points-1))); mean(w_dot_L_strk(2,1:(nr_points-1))); mean(w_dot_L_strk(3,1:(nr_points-1)))];
%     w_dot_R_strk_mean   = [mean(w_dot_R_strk(1,1:(nr_points-1))); mean(w_dot_R_strk(2,1:(nr_points-1))); mean(w_dot_R_strk(3,1:(nr_points-1)))];  
%     
%     for i = 1:nr_points
%         
%         w_dot_L_b(:,i)       = w_dot_L_b(:,i)-w_dot_L_b_mean;
%         w_dot_R_b(:,i)       = w_dot_R_b(:,i)-w_dot_R_b_mean;
%         w_dot_L_strk(:,i)    = w_dot_L_strk(:,i)-w_dot_L_strk_mean;
%         w_dot_R_strk(:,i)    = w_dot_R_strk(:,i)-w_dot_R_strk_mean;
%         
% %         wL_b(:,i)      = wL_b(:,i);
% %         wR_b(:,i)      = wR_b(:,i);
% %         wL_strk(:,i)   = wL_strk(:,i);
% %         wR_strk(:,i)   = wR_strk(:,i);
%     
%     end
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_L,'r',t,theta_R,'b',t,ones(nr_points,1)*mean(theta_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(theta_R(1:(nr_points))),'b')
%     title('wing kinematics')
%     subplot(3,1,2); plot(t,eta_L,'r',t,eta_R,'b',t,ones(nr_points,1)*mean(eta_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(eta_R(1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,phi_L,'r',t,phi_R,'b',t,ones(nr_points,1)*mean(phi_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(phi_R(1:(nr_points))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot_L,'r',t,theta_dot_R,'b',t,ones(nr_points,1)*mean(theta_dot_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(theta_dot_R(1:(nr_points))),'b')
%     title('euler angle derivatives')
%     subplot(3,1,2); plot(t,eta_dot_L,'r',t,eta_dot_R,'b',t,ones(nr_points,1)*mean(eta_dot_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(eta_dot_R(1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,phi_dot_L,'r',t,phi_dot_R,'b',t,ones(nr_points,1)*mean(phi_dot_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(phi_dot_R(1:(nr_points))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot_L,'r',t,theta_ddot_R,'b',t,ones(nr_points,1)*mean(theta_ddot_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(theta_ddot_R(1:(nr_points))),'b')
%     title('euler angle second derivative')
%     subplot(3,1,2); plot(t,eta_ddot_L,'r',t,eta_ddot_R,'b',t,ones(nr_points,1)*mean(eta_ddot_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(eta_ddot_R(1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,phi_ddot_L,'r',t,phi_ddot_R,'b',t,ones(nr_points,1)*mean(phi_ddot_L(1:(nr_points))),'r',t,ones(nr_points,1)*mean(phi_ddot_R(1:(nr_points))),'b')
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL(1,:),'r',t,wR(1,:),'b',t,ones(nr_points,1)*mean(wL(1,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR(1,1:(nr_points))),'b')
%     title('wL & wR')
%     subplot(3,1,2); plot(t,wL(2,:),'r',t,wR(2,:),'b',t,ones(nr_points,1)*mean(wL(2,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR(2,1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,wL(3,:),'r',t,wR(3,:),'b',t,ones(nr_points,1)*mean(wL(3,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR(3,1:(nr_points))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L(1,:),'r',t,w_dot_R(1,:),'b',t,ones(nr_points,1)*mean(w_dot_L(1,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R(1,1:(nr_points))),'b')
%     title('w dot L & w dot R')
%     subplot(3,1,2); plot(t,w_dot_L(2,:),'r',t,w_dot_R(2,:),'b',t,ones(nr_points,1)*mean(w_dot_L(2,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R(2,1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,w_dot_L(3,:),'r',t,w_dot_R(3,:),'b',t,ones(nr_points,1)*mean(w_dot_L(3,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R(3,1:(nr_points))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL_b(1,:),'r',t,wR_b(1,:),'b',t,ones(nr_points,1)*mean(wL_b(1,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR_b(1,1:(nr_points))),'b')
%     title('wL b & wR b')
%     subplot(3,1,2); plot(t,wL_b(2,:),'r',t,wR_b(2,:),'b',t,ones(nr_points,1)*mean(wL_b(2,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR_b(2,1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,wL_b(3,:),'r',t,wR_b(3,:),'b',t,ones(nr_points,1)*mean(wL_b(3,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR_b(3,1:(nr_points))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L_b(1,:),'r',t,w_dot_R_b(1,:),'b',t,ones(nr_points,1)*mean(w_dot_L_b(1,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R_b(1,1:(nr_points))),'b')
%     title('w dot L b & w dot R b')
%     subplot(3,1,2); plot(t,w_dot_L_b(2,:),'r',t,w_dot_R_b(2,:),'b',t,ones(nr_points,1)*mean(w_dot_L_b(2,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R_b(2,1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,w_dot_L_b(3,:),'r',t,w_dot_R_b(3,:),'b',t,ones(nr_points,1)*mean(w_dot_L_b(3,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R_b(3,1:(nr_points))),'b')
%     hold off    
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,wL_strk(1,:),'r',t,wR_strk(1,:),'b',t,ones(nr_points,1)*mean(wL_strk(1,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR_strk(1,1:(nr_points))),'b')
%     title('wL strk & wR strk')
%     subplot(3,1,2); plot(t,wL_strk(2,:),'r',t,wR_strk(2,:),'b',t,ones(nr_points,1)*mean(wL_strk(2,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR_strk(2,1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,wL_strk(3,:),'r',t,wR_strk(3,:),'b',t,ones(nr_points,1)*mean(wL_strk(3,1:(nr_points))),'r',t,ones(nr_points,1)*mean(wR_strk(3,1:(nr_points))),'b')
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,w_dot_L_strk(1,:),'r',t,w_dot_R_strk(1,:),'b',t,ones(nr_points,1)*mean(w_dot_L_strk(1,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R_strk(1,1:(nr_points))),'b')
%     title('w dot L strk & w dot R strk')
%     subplot(3,1,2); plot(t,w_dot_L_strk(2,:),'r',t,w_dot_R_strk(2,:),'b',t,ones(nr_points,1)*mean(w_dot_L_strk(2,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R_strk(2,1:(nr_points))),'b')
%     subplot(3,1,3); plot(t,w_dot_L_strk(3,:),'r',t,w_dot_R_strk(3,:),'b',t,ones(nr_points,1)*mean(w_dot_L_strk(3,1:(nr_points))),'r',t,ones(nr_points,1)*mean(w_dot_R_strk(3,1:(nr_points))),'b')
%     hold off    
    
%     Wt_L1 = [ 0.05; -1; 0];
%     Wt_L2 = [-0.05; -1 ; 0];
%     Wt_R1 = [ 0.05; 1; 0];
%     Wt_R2 = [-0.05; 1 ; 0];
%     
%     xyz_wing_L1 = zeros(3,nr_points);
%     xyz_wing_L2 = zeros(3,nr_points);
%     xyz_wing_R1 = zeros(3,nr_points);
%     xyz_wing_R2 = zeros(3,nr_points);
%     
%     for k = 1:nr_points
%         
%         xyz_wing_L1(:,k) = RL(:,:,k)'*Wt_L1;
%         xyz_wing_L2(:,k) = RL(:,:,k)'*Wt_L2;
%         xyz_wing_R1(:,k) = RR(:,:,k)'*Wt_R1;
%         xyz_wing_R2(:,k) = RR(:,:,k)'*Wt_R2;
%         
%     end
%     
%     alfa = 0:(2*pi/99):2*pi;
%     x_circ = cos(alfa);
%     y_circ = sin(alfa);
%     z_circ = zeros(1,100);
%     xyz_circ = [x_circ; y_circ; z_circ];
%     
%     xyz_strk = zeros(3,100);
%     
%     for l = 1:100
%         xyz_strk(:,l) = R_strk'*xyz_circ(:,l);
%     end
%     
%     figure()
%     hold on
%     for k = 1:nr_points
%         plot3([xyz_wing_L1(1,k) xyz_wing_L2(1,k)], [xyz_wing_L1(2,k) xyz_wing_L2(2,k)], [xyz_wing_L1(3,k) xyz_wing_L2(3,k)],'k')
%         plot3([xyz_wing_R1(1,k) xyz_wing_R2(1,k)], [xyz_wing_R1(2,k) xyz_wing_R2(2,k)], [xyz_wing_R1(3,k) xyz_wing_R2(3,k)],'k')
%         plot3(xyz_wing_L1(1,k),xyz_wing_L1(2,k),xyz_wing_L1(3,k),'o','Color','r')
%         plot3(xyz_wing_R1(1,k),xyz_wing_R1(2,k),xyz_wing_R1(3,k),'o','Color','g')
%     end
%     plot3(xyz_strk(1,:),xyz_strk(2,:),xyz_strk(3,:),'b')
%     plot3(xyz_wing_L1(1,1),xyz_wing_L1(2,1),xyz_wing_L1(3,1),'o','Color','b')
%     plot3([0 1],[0 0],[0 0],'r')
%     plot3([0 0],[0 0],[0 1],'b')
%     hold off
%     axis equal
%     pause

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
    kine.qL             = qL;
    kine.qR             = qR;
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



