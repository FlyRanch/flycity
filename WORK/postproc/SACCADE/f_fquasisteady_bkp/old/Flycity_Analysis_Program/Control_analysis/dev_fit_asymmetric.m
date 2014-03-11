function [ b, FM_LR, nr_wb_star ] = dev_fit_asymmetric( FM, a_dev_theta_L, a_dev_eta_L, a_dev_phi_L, a_dev_theta_R, a_dev_eta_R, a_dev_phi_R, down_up )

    % Find relations between a_dev and the maneuvering force or moment.
    
    % Initialization:
    
    n_pol_theta     = (size(a_dev_theta_L,1)-2)/2;
    n_pol_eta       = (size(a_dev_eta_L,1)-2)/2;
    n_pol_phi       = (size(a_dev_phi_L,1)-2)/2;

    % Sort FM:
    
    [FM_sort,FM_sort_id] = sort(FM);
    
    a_dev_theta_L_sort    = a_dev_theta_L(:,FM_sort_id);
    a_dev_eta_L_sort      = a_dev_eta_L(:,FM_sort_id);
    a_dev_phi_L_sort      = a_dev_phi_L(:,FM_sort_id);

    a_dev_theta_R_sort    = a_dev_theta_R(:,FM_sort_id);
    a_dev_eta_R_sort      = a_dev_eta_R(:,FM_sort_id);
    a_dev_phi_R_sort      = a_dev_phi_R(:,FM_sort_id);

    % Remove (weird) wingbeats with sum(a_dev) > 4*(standard deviation):
    
    [ FM_star, a_dev_theta_L_star, a_dev_eta_L_star, a_dev_phi_L_star, a_dev_theta_R_star, a_dev_eta_R_star, a_dev_phi_R_star, nr_wb_star ] = ...
    Remove_missed_wingbeats_asym( FM_sort, a_dev_theta_L_sort, a_dev_eta_L_sort, a_dev_phi_L_sort, a_dev_theta_R_sort, a_dev_eta_R_sort, a_dev_phi_R_sort );


    % Mirror left and right:

    split_id  = find(FM_star>=0,1,'first');

    if split_id == 1
        
        FM_LR = FM_star;

        a_dev_theta_L_t   = a_dev_theta_L_star;
        a_dev_eta_L_t     = a_dev_eta_L_star;
        a_dev_phi_L_t     = a_dev_phi_L_star;

        a_dev_theta_R_t   = a_dev_theta_R_star;
        a_dev_eta_R_t     = a_dev_eta_R_star;
        a_dev_phi_R_t     = a_dev_phi_R_star;

        b_theta_L_t     = lin_fit( a_dev_theta_L_t, FM_LR);
        b_eta_L_t       = lin_fit( a_dev_eta_L_t, FM_LR);
        b_phi_L_t       = lin_fit( a_dev_phi_L_t, FM_LR);

        b_theta_R_t     = lin_fit( a_dev_theta_R_t, FM_LR);
        b_eta_R_t       = lin_fit( a_dev_eta_R_t, FM_LR);
        b_phi_R_t       = lin_fit( a_dev_phi_R_t, FM_LR);

        [b_theta_L1,b_theta_L2] = average_fit(b_theta_L_t(1:(n_pol_theta+1),1),b_theta_L_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_L1,b_eta_L2] = average_fit(b_eta_L_t(1:(n_pol_eta+1),1),b_eta_L_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_L1,b_phi_L2] = average_fit(b_phi_L_t(1:(n_pol_phi+1),1),b_phi_L_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        [b_theta_R1,b_theta_R2] = average_fit(b_theta_R_t(1:(n_pol_theta+1),1),b_theta_R_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_R1,b_eta_R2] = average_fit(b_eta_R_t(1:(n_pol_eta+1),1),b_eta_R_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_R1,b_phi_R2] = average_fit(b_phi_R_t(1:(n_pol_phi+1),1),b_phi_R_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        b_theta_L   = [ [b_theta_L1;b_theta_L2] b_theta_L_t(:,2:3) ];
        b_eta_L     = [ [b_eta_L1;b_eta_L2] b_eta_L_t(:,2:3) ];
        b_phi_L     = [ [b_phi_L1;b_phi_L2] b_phi_L_t(:,2:3) ];

        b_theta_R   = [ [b_theta_R1;b_theta_R2] b_theta_R_t(:,2:3) ];
        b_eta_R     = [ [b_eta_R1;b_eta_R2] b_eta_R_t(:,2:3) ];
        b_phi_R     = [ [b_phi_R1;b_phi_R2] b_phi_R_t(:,2:3) ];
        
    elseif isempty(split_id) == 1
        
        FM_LR = -FM_star;

        a_dev_theta_L_t   = a_dev_theta_R_star;
        a_dev_eta_L_t     = a_dev_eta_R_star;
        a_dev_phi_L_t     = a_dev_phi_R_star;

        a_dev_theta_R_t   = a_dev_theta_L_star;
        a_dev_eta_R_t     = a_dev_eta_L_star;
        a_dev_phi_R_t     = a_dev_phi_L_star;

        b_theta_L_t     = lin_fit( a_dev_theta_L_t, FM_LR);
        b_eta_L_t       = lin_fit( a_dev_eta_L_t, FM_LR);
        b_phi_L_t       = lin_fit( a_dev_phi_L_t, FM_LR);

        b_theta_R_t     = lin_fit( a_dev_theta_R_t, FM_LR);
        b_eta_R_t       = lin_fit( a_dev_eta_R_t, FM_LR);
        b_phi_R_t       = lin_fit( a_dev_phi_R_t, FM_LR);

        [b_theta_L1,b_theta_L2] = average_fit(b_theta_L_t(1:(n_pol_theta+1),1),b_theta_L_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_L1,b_eta_L2] = average_fit(b_eta_L_t(1:(n_pol_eta+1),1),b_eta_L_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_L1,b_phi_L2] = average_fit(b_phi_L_t(1:(n_pol_phi+1),1),b_phi_L_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        [b_theta_R1,b_theta_R2] = average_fit(b_theta_R_t(1:(n_pol_theta+1),1),b_theta_R_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_R1,b_eta_R2] = average_fit(b_eta_R_t(1:(n_pol_eta+1),1),b_eta_R_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_R1,b_phi_R2] = average_fit(b_phi_R_t(1:(n_pol_phi+1),1),b_phi_R_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        b_theta_L   = [ [b_theta_L1;b_theta_L2] b_theta_L_t(:,2:3) ];
        b_eta_L     = [ [b_eta_L1;b_eta_L2] b_eta_L_t(:,2:3) ];
        b_phi_L     = [ [b_phi_L1;b_phi_L2] b_phi_L_t(:,2:3) ];

        b_theta_R   = [ [b_theta_R1;b_theta_R2] b_theta_R_t(:,2:3) ];
        b_eta_R     = [ [b_eta_R1;b_eta_R2] b_eta_R_t(:,2:3) ];
        b_phi_R     = [ [b_phi_R1;b_phi_R2] b_phi_R_t(:,2:3) ];
        
    else
        
        [FM_LR FM_LR_id] = sort([-FM_star(1:(split_id-1)) FM_star(split_id:nr_wb_star)]);

        a_dev_theta_L_t   = [a_dev_theta_R_star(:,1:(split_id-1)) a_dev_theta_L_star(:,split_id:nr_wb_star)];
        a_dev_eta_L_t     = [a_dev_eta_R_star(:,1:(split_id-1)) a_dev_eta_L_star(:,split_id:nr_wb_star)];
        a_dev_phi_L_t     = [a_dev_phi_R_star(:,1:(split_id-1)) a_dev_phi_L_star(:,split_id:nr_wb_star)];

        a_dev_theta_R_t   = [a_dev_theta_L_star(:,1:(split_id-1)) a_dev_theta_R_star(:,split_id:nr_wb_star)];
        a_dev_eta_R_t     = [a_dev_eta_L_star(:,1:(split_id-1)) a_dev_eta_R_star(:,split_id:nr_wb_star)];
        a_dev_phi_R_t     = [a_dev_phi_L_star(:,1:(split_id-1)) a_dev_phi_R_star(:,split_id:nr_wb_star)];

        b_theta_L_t     = lin_fit( a_dev_theta_L_t(:,FM_LR_id), FM_LR);
        b_eta_L_t       = lin_fit( a_dev_eta_L_t(:,FM_LR_id), FM_LR);
        b_phi_L_t       = lin_fit( a_dev_phi_L_t(:,FM_LR_id), FM_LR);

        b_theta_R_t     = lin_fit( a_dev_theta_R_t(:,FM_LR_id), FM_LR);
        b_eta_R_t       = lin_fit( a_dev_eta_R_t(:,FM_LR_id), FM_LR);
        b_phi_R_t       = lin_fit( a_dev_phi_R_t(:,FM_LR_id), FM_LR);

        [b_theta_L1,b_theta_L2] = average_fit(b_theta_L_t(1:(n_pol_theta+1),1),b_theta_L_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_L1,b_eta_L2] = average_fit(b_eta_L_t(1:(n_pol_eta+1),1),b_eta_L_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_L1,b_phi_L2] = average_fit(b_phi_L_t(1:(n_pol_phi+1),1),b_phi_L_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);

        [b_theta_R1,b_theta_R2] = average_fit(b_theta_R_t(1:(n_pol_theta+1),1),b_theta_R_t((n_pol_theta+2):(2*n_pol_theta+2),1),n_pol_theta,1,1,down_up);
        [b_eta_R1,b_eta_R2] = average_fit(b_eta_R_t(1:(n_pol_eta+1),1),b_eta_R_t((n_pol_eta+2):(2*n_pol_eta+2),1),n_pol_eta,1,1,down_up);
        [b_phi_R1,b_phi_R2] = average_fit(b_phi_R_t(1:(n_pol_phi+1),1),b_phi_R_t((n_pol_phi+2):(2*n_pol_phi+2),1),n_pol_phi,1,1,down_up);
        
        b_theta_L   = [ [b_theta_L1;b_theta_L2] b_theta_L_t(:,2:3) ];
        b_eta_L     = [ [b_eta_L1;b_eta_L2] b_eta_L_t(:,2:3) ];
        b_phi_L     = [ [b_phi_L1;b_phi_L2] b_phi_L_t(:,2:3) ];

        b_theta_R   = [ [b_theta_R1;b_theta_R2] b_theta_R_t(:,2:3) ];
        b_eta_R     = [ [b_eta_R1;b_eta_R2] b_eta_R_t(:,2:3) ];
        b_phi_R     = [ [b_phi_R1;b_phi_R2] b_phi_R_t(:,2:3) ];

    end

    b.theta_L   = b_theta_L;
    b.eta_L     = b_eta_L;
    b.phi_L     = b_phi_L;

    b.theta_R   = b_theta_R;
    b.eta_R     = b_eta_R;
    b.phi_R     = b_phi_R;    


end

