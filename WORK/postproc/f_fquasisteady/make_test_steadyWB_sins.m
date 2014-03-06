function [phi_L,phi_R,phi_dot_L,phi_dot_R,phi_ddot_L,phi_ddot_R,...
    theta_L,theta_R,theta_dot_L,theta_dot_R,theta_ddot_L,theta_ddot_R,...
    eta_L,eta_R,eta_dot_L,eta_dot_R,eta_ddot_L,eta_ddot_R,t,dt] = ...
    make_test_steadyWB_sins(nr_timepoints,f)


    theta_L0      = (5*pi)/180;
    eta_L0        = pi/2;
    phi_L0        = (10*pi)/180;
    theta_R0      = (5*pi)/180;
    eta_R0        = pi/2;
    phi_R0        = (10*pi)/180;
    A_theta_L     = (10*pi)/180;
    A_eta_L       = -(45*pi)/180;
    A_phi_L       = (60*pi)/180;
    A_theta_R     = (10*pi)/180;
    A_eta_R       = -(45*pi)/180;
    A_phi_R       = (60*pi)/180;

    
    
    dt = (1/f)/(nr_timepoints-1);
    t = 0:dt:(dt*(nr_timepoints-1));  

    theta_L         = (theta_L0+A_theta_L*cos(4*pi*f*t));
    eta_L           = -(eta_L0+A_eta_L*sin(2*pi*f*t));
    phi_L           = -(phi_L0+A_phi_L*cos(2*pi*f*t));
    
    theta_R         = -(theta_R0+A_theta_R*cos(4*pi*f*t));
    eta_R           = -(eta_R0+A_eta_R*sin(2*pi*f*t));
    phi_R           = (phi_R0+A_phi_R*cos(2*pi*f*t));

    theta_dot_L     = -A_theta_L*4*pi*f*sin(4*pi*f*t);
    eta_dot_L       = -A_eta_L*2*pi*f*cos(2*pi*f*t);
    phi_dot_L       = A_phi_L*2*pi*f*sin(2*pi*f*t);
    
    theta_dot_R     = A_theta_R*4*pi*f*sin(4*pi*f*t);
    eta_dot_R       = -A_eta_R*2*pi*f*cos(2*pi*f*t);
    phi_dot_R       = -A_phi_R*2*pi*f*sin(2*pi*f*t);
    
    theta_ddot_L    = -A_theta_L*16*pi^2*f^2*cos(4*pi*f*t);
    eta_ddot_L      = A_eta_L*4*pi^2*f^2*sin(2*pi*f*t);
    phi_ddot_L      = A_phi_L*4*pi^2*f^2*cos(2*pi*f*t);
    
    theta_ddot_R    = A_theta_R*16*pi^2*f^2*cos(4*pi*f*t);
    eta_ddot_R      = A_eta_R*4*pi^2*f^2*sin(2*pi*f*t);
    phi_ddot_R      = -A_phi_R*4*pi^2*f^2*cos(2*pi*f*t);

    
    
    