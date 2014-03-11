function [ FM_star, a_dev_theta_L_star, a_dev_eta_L_star, a_dev_phi_L_star, a_dev_theta_R_star, a_dev_eta_R_star, a_dev_phi_R_star, nr_wb_star ] = Remove_missed_wingbeats_asym( FM, a_dev_theta_L, a_dev_eta_L, a_dev_phi_L, a_dev_theta_R, a_dev_eta_R, a_dev_phi_R )

    nr_wb = length(FM);

    a_dev_theta_L_star    = [];
    a_dev_eta_L_star      = [];
    a_dev_phi_L_star      = [];
    a_dev_theta_R_star    = [];
    a_dev_eta_R_star      = [];
    a_dev_phi_R_star      = [];
    FM_star             = [];
    
    a_dev_theta_L_t   = a_dev_theta_L;
    a_dev_eta_L_t     = a_dev_eta_L;
    a_dev_phi_L_t     = a_dev_phi_L;
    a_dev_theta_R_t   = a_dev_theta_R;
    a_dev_eta_R_t     = a_dev_eta_R;
    a_dev_phi_R_t     = a_dev_phi_R;
    FM_t            = FM;
    
    nr_wb_t         = nr_wb;
    
    for k = 1:2
        
        if k < 2

            theta_L_sd_t    = std(sum(a_dev_theta_L_t));
            eta_L_sd_t      = std(sum(a_dev_eta_L_t));
            phi_L_sd_t      = std(sum(a_dev_phi_L_t));

            theta_R_sd_t    = std(sum(a_dev_theta_R_t));
            eta_R_sd_t      = std(sum(a_dev_eta_R_t));
            phi_R_sd_t      = std(sum(a_dev_phi_R_t));

            theta_L_id_weird = find( abs(sum(a_dev_theta_L_t)) > (4*theta_L_sd_t) );
            eta_L_id_weird   = find( abs(sum(a_dev_eta_L_t)) > (4*eta_L_sd_t) );
            phi_L_id_weird   = find( abs(sum(a_dev_phi_L_t)) > (4*phi_L_sd_t) );

            theta_R_id_weird = find( abs(sum(a_dev_theta_R_t)) > (4*theta_R_sd_t) );
            eta_R_id_weird   = find( abs(sum(a_dev_eta_R_t)) > (4*eta_R_sd_t) );
            phi_R_id_weird   = find( abs(sum(a_dev_phi_R_t)) > (4*phi_R_sd_t) );

            id_weird_t  = sort([theta_L_id_weird eta_L_id_weird phi_L_id_weird theta_R_id_weird eta_R_id_weird phi_R_id_weird]);
            id_weird    = unique(id_weird_t);

            a_dev_theta_L_t2  = [];
            a_dev_eta_L_t2    = [];
            a_dev_phi_L_t2    = [];
            a_dev_theta_R_t2  = [];
            a_dev_eta_R_t2    = [];
            a_dev_phi_R_t2    = [];
            FM_t2             = [];

            for i = 1:nr_wb_t

                if isempty(find(id_weird == i,1)) == 1

                    a_dev_theta_L_t2    = [a_dev_theta_L_t2 a_dev_theta_L_t(:,i)];
                    a_dev_eta_L_t2      = [a_dev_eta_L_t2 a_dev_eta_L_t(:,i)];
                    a_dev_phi_L_t2      = [a_dev_phi_L_t2 a_dev_phi_L_t(:,i)];
                    a_dev_theta_R_t2    = [a_dev_theta_R_t2 a_dev_theta_R_t(:,i)];
                    a_dev_eta_R_t2      = [a_dev_eta_R_t2 a_dev_eta_R_t(:,i)];
                    a_dev_phi_R_t2      = [a_dev_phi_R_t2 a_dev_phi_R_t(:,i)];
                    FM_t2               = [FM_t2 FM_t(i)];

                end

            end

            a_dev_theta_L_t   = a_dev_theta_L_t2;
            a_dev_eta_L_t     = a_dev_eta_L_t2;
            a_dev_phi_L_t     = a_dev_phi_L_t2;
            a_dev_theta_R_t   = a_dev_theta_R_t2;
            a_dev_eta_R_t     = a_dev_eta_R_t2;
            a_dev_phi_R_t     = a_dev_phi_R_t2;
            FM_t            = FM_t2;
            nr_wb_t         = length(FM_t2);
        
        else
            
            theta_L_sd_t    = std(sum(a_dev_theta_L_t));
            eta_L_sd_t      = std(sum(a_dev_eta_L_t));
            phi_L_sd_t      = std(sum(a_dev_phi_L_t));

            theta_R_sd_t    = std(sum(a_dev_theta_R_t));
            eta_R_sd_t      = std(sum(a_dev_eta_R_t));
            phi_R_sd_t      = std(sum(a_dev_phi_R_t));

            theta_L_id_weird = find( abs(sum(a_dev_theta_L_t)) > (4*theta_L_sd_t) );
            eta_L_id_weird   = find( abs(sum(a_dev_eta_L_t)) > (4*eta_L_sd_t) );
            phi_L_id_weird   = find( abs(sum(a_dev_phi_L_t)) > (4*phi_L_sd_t) );

            theta_R_id_weird = find( abs(sum(a_dev_theta_R_t)) > (4*theta_R_sd_t) );
            eta_R_id_weird   = find( abs(sum(a_dev_eta_R_t)) > (4*eta_R_sd_t) );
            phi_R_id_weird   = find( abs(sum(a_dev_phi_R_t)) > (4*phi_R_sd_t) );

            id_weird_t  = sort([theta_L_id_weird eta_L_id_weird phi_L_id_weird theta_R_id_weird eta_R_id_weird phi_R_id_weird]);
            id_weird    = unique(id_weird_t);

            a_dev_theta_L_t2  = [];
            a_dev_eta_L_t2    = [];
            a_dev_phi_L_t2    = [];
            a_dev_theta_R_t2  = [];
            a_dev_eta_R_t2    = [];
            a_dev_phi_R_t2    = [];
            FM_t2             = [];

            for i = 1:nr_wb_t

                if isempty(find(id_weird == i,1)) == 1

                    a_dev_theta_L_t2    = [a_dev_theta_L_t2 a_dev_theta_L_t(:,i)];
                    a_dev_eta_L_t2      = [a_dev_eta_L_t2 a_dev_eta_L_t(:,i)];
                    a_dev_phi_L_t2      = [a_dev_phi_L_t2 a_dev_phi_L_t(:,i)];
                    a_dev_theta_R_t2    = [a_dev_theta_R_t2 a_dev_theta_R_t(:,i)];
                    a_dev_eta_R_t2      = [a_dev_eta_R_t2 a_dev_eta_R_t(:,i)];
                    a_dev_phi_R_t2      = [a_dev_phi_R_t2 a_dev_phi_R_t(:,i)];
                    FM_t2               = [FM_t2 FM_t(i)];

                end

            end

            a_dev_theta_L_star   = a_dev_theta_L_t2;
            a_dev_eta_L_star     = a_dev_eta_L_t2;
            a_dev_phi_L_star     = a_dev_phi_L_t2;
            a_dev_theta_R_star   = a_dev_theta_R_t2;
            a_dev_eta_R_star     = a_dev_eta_R_t2;
            a_dev_phi_R_star     = a_dev_phi_R_t2;
            FM_star              = FM_t2;
            nr_wb_star           = length(FM_t2);      
            
        end
        
    end
    

end

