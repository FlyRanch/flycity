function [ FM_star, a_dev_theta_star, a_dev_eta_star, a_dev_phi_star, nr_wb_star ] = Remove_missed_wingbeats_sym( FM, a_dev_theta, a_dev_eta, a_dev_phi )

    nr_wb = length(FM);
    
%     theta_mean  = mean(sum(a_dev_theta));
%     eta_mean    = mean(sum(a_dev_eta));
%     phi_mean    = mean(sum(a_dev_phi));
%     
%     theta_sd    = std(sum(a_dev_theta));
%     eta_sd      = std(sum(a_dev_eta));
%     phi_sd      = std(sum(a_dev_phi));
%     
%     figure()
%     plot(1:nr_wb,FM)
%     
%     figure()
%     plot(1:nr_wb,sum(a_dev_theta),1:nr_wb,ones(nr_wb,1)*theta_mean,1:nr_wb,ones(nr_wb,1)*(theta_mean+theta_sd),1:nr_wb,ones(nr_wb,1)*(theta_mean-theta_sd))
%     
%     figure()
%     plot(1:nr_wb,sum(a_dev_eta),1:nr_wb,ones(nr_wb,1)*eta_mean,1:nr_wb,ones(nr_wb,1)*(eta_mean+eta_sd),1:nr_wb,ones(nr_wb,1)*(eta_mean-eta_sd))
%     
%     figure()
%     plot(1:nr_wb,sum(a_dev_phi),1:nr_wb,ones(nr_wb,1)*phi_mean,1:nr_wb,ones(nr_wb,1)*(phi_mean+phi_sd),1:nr_wb,ones(nr_wb,1)*(phi_mean-phi_sd))
%     
    a_dev_theta_star    = [];
    a_dev_eta_star      = [];
    a_dev_phi_star      = [];
    FM_star             = [];
    
    a_dev_theta_t   = a_dev_theta;
    a_dev_eta_t     = a_dev_eta;
    a_dev_phi_t     = a_dev_phi;
    FM_t            = FM;
    
    nr_wb_t         = nr_wb;
    
    for k = 1:2
        
        if k < 2

            theta_sd_t    = std(sum(a_dev_theta_t));
            eta_sd_t      = std(sum(a_dev_eta_t));
            phi_sd_t      = std(sum(a_dev_phi_t));

            theta_id_weird = find( abs(sum(a_dev_theta_t)) > (4*theta_sd_t) );
            eta_id_weird   = find( abs(sum(a_dev_eta_t)) > (4*eta_sd_t) );
            phi_id_weird   = find( abs(sum(a_dev_phi_t)) > (4*phi_sd_t) );

            id_weird_t  = sort([theta_id_weird eta_id_weird phi_id_weird]);
            id_weird    = unique(id_weird_t);

            a_dev_theta_t2  = [];
            a_dev_eta_t2    = [];
            a_dev_phi_t2    = [];
            FM_t2           = [];

            for i = 1:nr_wb_t

                if isempty(find(id_weird == i,1)) == 1

                    a_dev_theta_t2    = [a_dev_theta_t2 a_dev_theta_t(:,i)];
                    a_dev_eta_t2      = [a_dev_eta_t2 a_dev_eta_t(:,i)];
                    a_dev_phi_t2      = [a_dev_phi_t2 a_dev_phi_t(:,i)];
                    FM_t2             = [FM_t2 FM_t(i)];

                end

            end

            a_dev_theta_t   = a_dev_theta_t2;
            a_dev_eta_t     = a_dev_eta_t2;
            a_dev_phi_t     = a_dev_phi_t2;
            FM_t            = FM_t2;
            nr_wb_t         = length(FM_t2);
        
        else
            
            theta_sd_t    = std(sum(a_dev_theta_t));
            eta_sd_t      = std(sum(a_dev_eta_t));
            phi_sd_t      = std(sum(a_dev_phi_t));

            theta_id_weird = find( abs(sum(a_dev_theta_t)) > (4*theta_sd_t) );
            eta_id_weird   = find( abs(sum(a_dev_eta_t)) > (4*eta_sd_t) );
            phi_id_weird   = find( abs(sum(a_dev_phi_t)) > (4*phi_sd_t) );

            id_weird_t  = sort([theta_id_weird eta_id_weird phi_id_weird]);
            id_weird    = unique(id_weird_t);

            a_dev_theta_t2  = [];
            a_dev_eta_t2    = [];
            a_dev_phi_t2    = [];
            FM_t2           = [];

            for i = 1:nr_wb_t

                if isempty(find(id_weird == i,1)) == 1

                    a_dev_theta_t2    = [a_dev_theta_t2 a_dev_theta_t(:,i)];
                    a_dev_eta_t2      = [a_dev_eta_t2 a_dev_eta_t(:,i)];
                    a_dev_phi_t2      = [a_dev_phi_t2 a_dev_phi_t(:,i)];
                    FM_t2             = [FM_t2 FM_t(i)];

                end

            end

            a_dev_theta_star   = a_dev_theta_t2;
            a_dev_eta_star     = a_dev_eta_t2;
            a_dev_phi_star     = a_dev_phi_t2;
            FM_star            = FM_t2;
            nr_wb_star         = length(FM_t2);            
            
        end
        
    end
    
%     figure()
%     plot(1:nr_wb_star,FM_star)
%     
%     figure()
%     plot(1:nr_wb_star,sum(a_dev_theta_star))
%     
%     figure()
%     plot(1:nr_wb_star,sum(a_dev_eta_star))
%     
%     figure()
%     plot(1:nr_wb_star,sum(a_dev_phi_star))
    

end

