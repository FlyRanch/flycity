function [ b, FM_n, FM_p, N_n, N_p, N ] = dev_fit_symmetric( FM, a_dev_theta, a_dev_eta, a_dev_phi, down_up )

    % Find relations between a_dev and the maneuvering force or moment.
    
    % Initialization:
    
    nr_wb = size(a_dev_theta,2);
    
    n_pol_theta     = (size(a_dev_theta,1)-2)/2;
    n_pol_eta       = (size(a_dev_eta,1)-2)/2;
    n_pol_phi       = (size(a_dev_phi,1)-2)/2;
    
    % Sort FM:
    
    [ FM_sort, FM_sort_id ] = sort(FM);
    
    a_dev_theta_sort    = a_dev_theta(:,FM_sort_id);
    a_dev_eta_sort      = a_dev_eta(:,FM_sort_id);
    a_dev_phi_sort      = a_dev_phi(:,FM_sort_id);
    
    
    % Remove (weird) wingbeats with sum(a_dev) > 4*(standard deviation):
    
    [ FM_star, a_dev_theta_star, a_dev_eta_star, a_dev_phi_star, nr_wb_star ] = Remove_missed_wingbeats_sym( FM_sort, a_dev_theta_sort, a_dev_eta_sort, a_dev_phi_sort );


    split_id  = find(FM_star>=0,1,'first');
    
    N_n     = 0;
    N_p     = 0;
    N       = 0;

    if split_id == 1

        FM_p = FM_star;
        FM_n = [];

        a_dev_theta_p_t = a_dev_theta_star;
        a_dev_eta_p_t   = a_dev_eta_star;
        a_dev_phi_p_t   = a_dev_phi_star;

        b_theta_p_t     = lin_fit( a_dev_theta_p_t, FM_p );
        b_eta_p_t       = lin_fit( a_dev_eta_p_t, FM_p );
        b_phi_p_t       = lin_fit( a_dev_phi_p_t, FM_p );

        [b_theta_p1,b_theta_p2] = average_fit(b_theta_p_t(1:(n_pol_theta+1),1,k),b_theta_p_t((n_pol_theta+2):(2*n_pol_theta+2),1,k),n_pol_theta,1,1,down_up);
        [b_eta_p1,b_eta_p2] = average_fit(b_eta_p_t(1:(n_pol_eta+1),1,k),b_eta_p_t((n_pol_eta+2):(2*n_pol_eta+2),1,k),n_pol_eta,1,1,down_up);
        [b_phi_p1,b_phi_p2] = average_fit(b_phi_p_t(1:(n_pol_phi+1),1,k),b_phi_p_t((n_pol_phi+2):(2*n_pol_phi+2),1,k),n_pol_phi,1,1,down_up);

        b_theta_p   = [ [b_theta_p1;b_theta_p2] b_theta_p_t(:,2:3) ];
        b_eta_p     = [ [b_eta_p1;b_eta_p2] b_eta_p_t(:,2:3) ];
        b_phi_p     = [ [b_phi_p1;b_phi_p2] b_phi_p_t(:,2:3) ];

        b_theta_n     = zeros(n_pol_theta*2+2,3);
        b_eta_n       = zeros(n_pol_eta*2+2,3);
        b_phi_n       = zeros(n_pol_phi*2+2,3);
        
        N_p     = nr_wb_star;
        N       = nr_wb_star;

    elseif isempty(split_id) == 1

        FM_p = [];
        FM_n = FM_star;

        a_dev_theta_n_t = a_dev_theta_star;
        a_dev_eta_n_t   = a_dev_eta_star;
        a_dev_phi_n_t   = a_dev_phi_star;

        b_theta_n_t     = lin_fit( a_dev_theta_n_t, FM_n );
        b_eta_n_t       = lin_fit( a_dev_eta_n_t, FM_n );
        b_phi_n_t       = lin_fit( a_dev_phi_n_t, FM_n );

        [b_theta_n1,b_theta_n2] = average_fit(b_theta_n_t(1:(n_pol_theta+1),1,k),b_theta_n_t((n_pol_theta+2):(2*n_pol_theta+2),1,k),n_pol_theta,1,1,down_up);
        [b_eta_n1,b_eta_n2] = average_fit(b_eta_n_t(1:(n_pol_eta+1),1,k),b_eta_n_t((n_pol_eta+2):(2*n_pol_eta+2),1,k),n_pol_eta,1,1,down_up);
        [b_phi_n1,b_phi_n2] = average_fit(b_phi_n_t(1:(n_pol_phi+1),1,k),b_phi_n_t((n_pol_phi+2):(2*n_pol_phi+2),1,k),n_pol_phi,1,1,down_up);

        b_theta_p   = zeros(n_pol_theta*2+2,3);
        b_eta_p     = zeros(n_pol_eta*2+2,3);
        b_phi_p     = zeros(n_pol_phi*2+2,3);

        b_theta_n   = [ [b_theta_n1;b_theta_n2] b_theta_n_t(:,2:3) ];
        b_eta_n     = [ [b_eta_n1;b_eta_n2] b_eta_n_t(:,2:3) ];
        b_phi_n     = [ [b_phi_n1;b_phi_n2] b_phi_n_t(:,2:3) ];

        N_n     = nr_wb_star;
        N       = nr_wb_star;

    else
        
        FM_p = FM_star(split_id:nr_wb_star);
        FM_n = FM_star(1:(split_id-1));

        a_dev_theta_p_t = a_dev_theta_star(:,split_id:nr_wb_star);
        a_dev_eta_p_t   = a_dev_eta_star(:,split_id:nr_wb_star);
        a_dev_phi_p_t   = a_dev_phi_star(:,split_id:nr_wb_star);

        a_dev_theta_n_t = a_dev_theta_star(:,1:(split_id-1));
        a_dev_eta_n_t   = a_dev_eta_star(:,1:(split_id-1));
        a_dev_phi_n_t   = a_dev_phi_star(:,1:(split_id-1));

        b_theta_p_t     = lin_fit( a_dev_theta_p_t, FM_p );
        b_eta_p_t       = lin_fit( a_dev_eta_p_t, FM_p );
        b_phi_p_t       = lin_fit( a_dev_phi_p_t, FM_p );

        b_theta_n_t     = lin_fit( a_dev_theta_n_t, FM_n );
        b_eta_n_t       = lin_fit( a_dev_eta_n_t, FM_n );
        b_phi_n_t       = lin_fit( a_dev_phi_n_t, FM_n );
        
        [b_theta_p1,b_theta_p2] = average_fit(b_theta_p_t(1:(n_pol_theta+1),1),b_theta_p_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_p1,b_eta_p2] = average_fit(b_eta_p_t(1:(n_pol_eta+1),1),b_eta_p_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_p1,b_phi_p2] = average_fit(b_phi_p_t(1:(n_pol_phi+1),1),b_phi_p_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        [b_theta_n1,b_theta_n2] = average_fit(b_theta_n_t(1:(n_pol_theta+1),1),b_theta_n_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_n1,b_eta_n2] = average_fit(b_eta_n_t(1:(n_pol_eta+1),1),b_eta_n_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_n1,b_phi_n2] = average_fit(b_phi_n_t(1:(n_pol_phi+1),1),b_phi_n_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        b_theta_p   = [ [b_theta_p1;b_theta_p2] b_theta_p_t(:,2:3) ];
        b_eta_p     = [ [b_eta_p1;b_eta_p2] b_eta_p_t(:,2:3) ];
        b_phi_p     = [ [b_phi_p1;b_phi_p2] b_phi_p_t(:,2:3) ];

        b_theta_n   = [ [b_theta_n1;b_theta_n2] b_theta_n_t(:,2:3) ];
        b_eta_n     = [ [b_eta_n1;b_eta_n2] b_eta_n_t(:,2:3) ];
        b_phi_n     = [ [b_phi_n1;b_phi_n2] b_phi_n_t(:,2:3) ];

        N_p     = split_id-1;
        N_n     = nr_wb_star-split_id+1;
        N       = nr_wb_star;

    end
    
    b.theta_p   = b_theta_p;
    b.eta_p     = b_eta_p;
    b.phi_p     = b_phi_p;

    b.theta_n   = b_theta_n;
    b.eta_n     = b_eta_n;
    b.phi_n     = b_phi_n;

end

