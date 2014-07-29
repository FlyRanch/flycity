function [ kine ] = poly_derivatives( a_fit, nr_points )


    % Use the polynomial fitting coefficients to compute: -----------------
    %
    %   theta, eta, phi
    %   theta_dot, eta_dot, phi_dot
    %   theta_ddot, eta_ddot, phi_ddot
    %
    % ---------------------------------------------------------------------

    
    a_theta         = a_fit.a_theta;
    a_eta           = a_fit.a_eta;
    a_phi           = a_fit.a_phi;
    
    f               = a_fit.f;
    
    down_up         = a_fit.down_up;
    
    n_pol_theta     = (length(a_theta)-2)/2;
    n_pol_eta       = (length(a_eta)-2)/2;
    n_pol_phi       = (length(a_phi)-2)/2;
    
    dt = (1/f)/(nr_points-1);
    
    nr_points_down  = round(down_up*nr_points);
    
    nr_points_up    = nr_points-nr_points_down;
    
    [ t, X_theta ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up, nr_points, 0, dt*(nr_points-1), 0 );
    [ ~, X_eta ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up, nr_points, 0, dt*(nr_points-1), 0 );
    [ ~, X_phi ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up, nr_points, 0, dt*(nr_points-1), 0 );
    
    [ ~, X_theta_dot ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up, nr_points, 0, dt*(nr_points-1), 1 );
    [ ~, X_eta_dot ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up, nr_points, 0, dt*(nr_points-1), 1 );
    [ ~, X_phi_dot ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up, nr_points, 0, dt*(nr_points-1), 1 );
    
    [ ~, X_theta_ddot ]  = Wingbeat_Legendre_matrix( n_pol_theta, down_up, nr_points, 0, dt*(nr_points-1), 2 );
    [ ~, X_eta_ddot ]    = Wingbeat_Legendre_matrix( n_pol_eta, down_up, nr_points, 0, dt*(nr_points-1), 2 );
    [ ~, X_phi_ddot ]    = Wingbeat_Legendre_matrix( n_pol_phi, down_up, nr_points, 0, dt*(nr_points-1), 2 );
    
    theta           = X_theta*a_theta;
    eta             = X_eta*a_eta;
    phi             = X_phi*a_phi;
    
    theta_dot_t       = (X_theta_dot*a_theta);
    eta_dot_t         = (X_eta_dot*a_eta);
    phi_dot_t         = (X_phi_dot*a_phi);
    
    theta_ddot_t       = (X_theta_ddot*a_theta);
    eta_ddot_t         = (X_eta_ddot*a_eta);
    phi_ddot_t         = (X_phi_ddot*a_phi);
    
    d_theta_dt = gradient(theta)./dt;
    d_eta_dt   = gradient(eta)./dt;
    d_phi_dt   = gradient(phi)./dt;
    
%     theta_dot(1:nr_points_down)     = theta_dot_t(1:nr_points_down).*(2*f/down_up);
%     eta_dot(1:nr_points_down)       = eta_dot_t(1:nr_points_down).*(2*f/down_up);
%     phi_dot(1:nr_points_down)       = phi_dot_t(1:nr_points_down).*(2*f/down_up);
%        
%     theta_dot((nr_points_down+1):nr_points)     = theta_dot_t((nr_points_down+1):nr_points).*(2*f/(1-down_up));
%     eta_dot((nr_points_down+1):nr_points)       = eta_dot_t((nr_points_down+1):nr_points).*(2*f/(1-down_up));
%     phi_dot((nr_points_down+1):nr_points)       = phi_dot_t((nr_points_down+1):nr_points).*(2*f/(1-down_up));
%     
%     theta_dot_diff  = theta_dot(nr_points_down)-theta_dot(nr_points_down+1)
%     eta_dot_diff    = eta_dot(nr_points_down)-eta_dot(nr_points_down+1)
%     phi_dot_diff    = phi_dot(nr_points_down)-phi_dot(nr_points_down+1)
%     
% %     theta_dot(1:nr_points_down)     = theta_dot(1:nr_points_down)-0.5*theta_dot_diff;
% %     eta_dot(1:nr_points_down)       = eta_dot(1:nr_points_down)-0.5*eta_dot_diff;
% %     phi_dot(1:nr_points_down)       = phi_dot(1:nr_points_down)-0.5*phi_dot_diff;
%     
%     theta_dot((nr_points_down+1):nr_points)     = theta_dot((nr_points_down+1):nr_points)+theta_dot_diff;
%     eta_dot((nr_points_down+1):nr_points)       = eta_dot((nr_points_down+1):nr_points)+eta_dot_diff;
%     phi_dot((nr_points_down+1):nr_points)       = phi_dot((nr_points_down+1):nr_points)+phi_dot_diff;
%        
%     
%     theta_ddot(1:nr_points_down)     = theta_ddot_t(1:nr_points_down).*(2*f/down_up);
%     eta_ddot(1:nr_points_down)       = eta_ddot_t(1:nr_points_down).*(2*f/down_up);
%     phi_ddot(1:nr_points_down)       = phi_ddot_t(1:nr_points_down).*(2*f/down_up);
%     
%     theta_ddot((nr_points_down+1):nr_points)     = theta_ddot_t((nr_points_down+1):nr_points).*(2*f/(1-down_up));
%     eta_ddot((nr_points_down+1):nr_points)       = eta_ddot_t((nr_points_down+1):nr_points).*(2*f/(1-down_up));
%     phi_ddot((nr_points_down+1):nr_points)       = phi_ddot_t((nr_points_down+1):nr_points).*(2*f/(1-down_up));

   

%     theta_dot(1:nr_points_down)     = theta_dot_t(1:nr_points_down).*(4*f);
%     eta_dot(1:nr_points_down)       = eta_dot_t(1:nr_points_down).*(4*f);
%     phi_dot(1:nr_points_down)       = phi_dot_t(1:nr_points_down).*(4*f);
%     
%     theta_dot((nr_points_down+1):nr_points)     = theta_dot_t((nr_points_down+1):nr_points).*(4*f);
%     eta_dot((nr_points_down+1):nr_points)       = eta_dot_t((nr_points_down+1):nr_points).*(4*f);
%     phi_dot((nr_points_down+1):nr_points)       = phi_dot_t((nr_points_down+1):nr_points).*(4*f);
%     
%     theta_ddot(1:nr_points_down)     = theta_ddot_t(1:nr_points_down).*(4*f);
%     eta_ddot(1:nr_points_down)       = eta_ddot_t(1:nr_points_down).*(4*f);
%     phi_ddot(1:nr_points_down)       = phi_ddot_t(1:nr_points_down).*(4*f);
%     
%     theta_ddot((nr_points_down+1):nr_points)     = theta_ddot_t((nr_points_down+1):nr_points).*(4*f);
%     eta_ddot((nr_points_down+1):nr_points)       = eta_ddot_t((nr_points_down+1):nr_points).*(4*f);
%     phi_ddot((nr_points_down+1):nr_points)       = phi_ddot_t((nr_points_down+1):nr_points).*(4*f);

    
    % Scale matrix:
    
    scale = zeros(nr_points,1);
    
    scale_amp_down  = ((0.5/down_up)-1);
    scale_amp_up    = ((0.5/(1-down_up))-1);
    
    scale(1:nr_points_down) = 1+scale_amp_down*sin(0:(pi/(nr_points_down-1)):pi).^0.25;
    scale(nr_points_down:nr_points) = 1+scale_amp_up*sin(0:(pi/(nr_points_up)):pi).^0.25;
    
%     scale(1:nr_points_down) = 1+scale_amp_down*sin(0:(pi/(nr_points_down-1)):pi).^0.01;
%     scale(nr_points_down:nr_points) = 1+scale_amp_up*sin(pi:(pi/(nr_points_up)):(2*pi)).^0.01;
   
%     figure()
%     plot(scale)
%     
%     mean(scale)
%     
%     pause
    
    theta_dot(1:nr_points_down)     = theta_dot_t(1:nr_points_down).*(4*f*scale(1:nr_points_down));
    eta_dot(1:nr_points_down)       = eta_dot_t(1:nr_points_down).*(4*f*scale(1:nr_points_down));
    phi_dot(1:nr_points_down)       = phi_dot_t(1:nr_points_down).*(4*f*scale(1:nr_points_down));
    
    theta_dot(nr_points_down:nr_points)     = theta_dot_t(nr_points_down:nr_points).*(4*f*scale(nr_points_down:nr_points));
    eta_dot(nr_points_down:nr_points)       = eta_dot_t(nr_points_down:nr_points).*(4*f*scale(nr_points_down:nr_points));
    phi_dot(nr_points_down:nr_points)       = phi_dot_t(nr_points_down:nr_points).*(4*f*scale(nr_points_down:nr_points));
    
    theta_ddot(1:nr_points_down)     = theta_ddot_t(1:nr_points_down).*(4*f*scale(1:nr_points_down));
    eta_ddot(1:nr_points_down)       = eta_ddot_t(1:nr_points_down).*(4*f*scale(1:nr_points_down));
    phi_ddot(1:nr_points_down)       = phi_ddot_t(1:nr_points_down).*(4*f*scale(1:nr_points_down));
    
    theta_ddot(nr_points_down:nr_points)     = theta_ddot_t(nr_points_down:nr_points).*(4*f*scale(nr_points_down:nr_points));
    eta_ddot(nr_points_down:nr_points)       = eta_ddot_t(nr_points_down:nr_points).*(4*f*scale(nr_points_down:nr_points));
    phi_ddot(nr_points_down:nr_points)       = phi_ddot_t(nr_points_down:nr_points).*(4*f*scale(nr_points_down:nr_points));
    
    % Compute wingbeat averages of derivatives and set to zero:
    
    theta_dot_mean      = mean(theta_dot(1:(nr_points-1)));
    eta_dot_mean        = mean(eta_dot(1:(nr_points-1)));
    phi_dot_mean        = mean(phi_dot(1:(nr_points-1)));
    
    theta_ddot_mean     = mean(theta_ddot(1:(nr_points-1)));
    eta_ddot_mean       = mean(eta_ddot(1:(nr_points-1)));
    phi_ddot_mean       = mean(phi_ddot(1:(nr_points-1)));
    
%     theta_dot  = theta_dot-theta_dot_mean;
%     eta_dot    = eta_dot-eta_dot_mean;
%     phi_dot    = phi_dot-phi_dot_mean;
%     
%     theta_ddot = theta_ddot-theta_ddot_mean;
%     eta_ddot   = eta_ddot-eta_ddot_mean;
%     phi_ddot   = phi_ddot-phi_ddot_mean;
    
    kine.t      = t;
    
    kine.theta  = theta;
    kine.eta    = eta;
    kine.phi    = phi;
    
    kine.theta_dot  = theta_dot;
    kine.eta_dot    = eta_dot;
    kine.phi_dot    = phi_dot;
    
    kine.theta_ddot = theta_ddot;
    kine.eta_ddot   = eta_ddot;
    kine.phi_ddot   = phi_ddot;
    

    
%     figure()
%     plot(wing_kin_norm)
% %     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta,t,ones(nr_points,1)*mean(theta))
%     subplot(3,1,2); plot(t,eta,t,ones(nr_points,1)*mean(eta))
%     subplot(3,1,3); plot(t,phi,t,ones(nr_points,1)*mean(phi))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_dot,'b',t,ones(nr_points,1)*mean(theta_dot(1:(nr_points-1))),'k',t,d_theta_dt,'r',t,ones(nr_points,1)*mean(d_theta_dt(1:(nr_points-1))),'g')
%     subplot(3,1,2); plot(t,eta_dot,'b',t,ones(nr_points,1)*mean(eta_dot(1:(nr_points-1))),'k',t,d_eta_dt,'r',t,ones(nr_points,1)*mean(d_eta_dt(1:(nr_points-1))),'g')
%     subplot(3,1,3); plot(t,phi_dot,'b',t,ones(nr_points,1)*mean(phi_dot(1:(nr_points-1))),'k',t,d_phi_dt,'r',t,ones(nr_points,1)*mean(d_phi_dt(1:(nr_points-1))),'g')
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(t,theta_ddot,t,ones(nr_points,1)*mean(theta_ddot(1:(nr_points-1))))
%     subplot(3,1,2); plot(t,eta_ddot,t,ones(nr_points,1)*mean(eta_ddot(1:(nr_points-1))))
%     subplot(3,1,3); plot(t,phi_ddot,t,ones(nr_points,1)*mean(phi_ddot(1:(nr_points-1))),t(nr_points_down),phi_ddot(nr_points_down),'o')
%     hold off
%     
%     pause

end

