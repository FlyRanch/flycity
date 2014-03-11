function [ weird_wb ] = Remove_missed_wingbeats( a_dev_theta_L, a_dev_eta_L, a_dev_phi_L, a_dev_theta_R, a_dev_eta_R, a_dev_phi_R )

    nr_wb = size(a_dev_theta_L,2);

%     a_dev_theta_L_star      = [];
%     a_dev_eta_L_star        = [];
%     a_dev_phi_L_star        = [];
%     a_dev_theta_R_star      = [];
%     a_dev_eta_R_star        = [];
%     a_dev_phi_R_star        = [];
    
    a_dev_theta_L_t         = a_dev_theta_L;
    a_dev_eta_L_t           = a_dev_eta_L;
    a_dev_phi_L_t           = a_dev_phi_L;
    a_dev_theta_R_t         = a_dev_theta_R;
    a_dev_eta_R_t           = a_dev_eta_R;
    a_dev_phi_R_t           = a_dev_phi_R;
    
    n_pol_theta             = (size(a_dev_theta_L,1)-2)/2;
    n_pol_eta               = (size(a_dev_eta_L,1)-2)/2;
    n_pol_phi               = (size(a_dev_phi_L,1)-2)/2;
    
    nr_wb_t                 = nr_wb;
    
    weird_wb                = [];
    
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
            
            weird_wb = id_weird;
            
            weird_wb
            
            a_dev_theta_L_t2  = [];
            a_dev_eta_L_t2    = [];
            a_dev_phi_L_t2    = [];
            a_dev_theta_R_t2  = [];
            a_dev_eta_R_t2    = [];
            a_dev_phi_R_t2    = [];

            for i = 1:nr_wb_t

                if isempty(find(id_weird == i,1)) == 0

                    a_dev_theta_L_t(:,i)    = nan(n_pol_theta*2+2,1);
                    a_dev_eta_L_t(:,i)      = nan(n_pol_eta*2+2,1);
                    a_dev_phi_L_t(:,i)      = nan(n_pol_phi*2+2,1);
                    a_dev_theta_R_t(:,i)    = nan(n_pol_theta*2+2,1);
                    a_dev_eta_R_t(:,i)      = nan(n_pol_eta*2+2,1);
                    a_dev_phi_R_t(:,i)      = nan(n_pol_phi*2+2,1);

                end
                
                if isempty(find(id_weird == i,1)) == 1

                    a_dev_theta_L_t2    = [a_dev_theta_L_t2 a_dev_theta_L_t(:,i)];
                    a_dev_eta_L_t2      = [a_dev_eta_L_t2 a_dev_eta_L_t(:,i)];
                    a_dev_phi_L_t2      = [a_dev_phi_L_t2 a_dev_phi_L_t(:,i)];
                    a_dev_theta_R_t2    = [a_dev_theta_R_t2 a_dev_theta_R_t(:,i)];
                    a_dev_eta_R_t2      = [a_dev_eta_R_t2 a_dev_eta_R_t(:,i)];
                    a_dev_phi_R_t2      = [a_dev_phi_R_t2 a_dev_phi_R_t(:,i)];

                end

            end
        
        else
            
            theta_L_sd_t    = nanstd(sum(a_dev_theta_L_t2));
            eta_L_sd_t      = nanstd(sum(a_dev_eta_L_t2));
            phi_L_sd_t      = nanstd(sum(a_dev_phi_L_t2));

            theta_R_sd_t    = nanstd(sum(a_dev_theta_R_t2));
            eta_R_sd_t      = nanstd(sum(a_dev_eta_R_t2));
            phi_R_sd_t      = nanstd(sum(a_dev_phi_R_t2));

            theta_L_id_weird = find( abs(sum(a_dev_theta_L_t)) > (4*theta_L_sd_t) );
            eta_L_id_weird   = find( abs(sum(a_dev_eta_L_t)) > (4*eta_L_sd_t) );
            phi_L_id_weird   = find( abs(sum(a_dev_phi_L_t)) > (4*phi_L_sd_t) );

            theta_R_id_weird = find( abs(sum(a_dev_theta_R_t)) > (4*theta_R_sd_t) );
            eta_R_id_weird   = find( abs(sum(a_dev_eta_R_t)) > (4*eta_R_sd_t) );
            phi_R_id_weird   = find( abs(sum(a_dev_phi_R_t)) > (4*phi_R_sd_t) );

            id_weird_t  = sort([theta_L_id_weird eta_L_id_weird phi_L_id_weird theta_R_id_weird eta_R_id_weird phi_R_id_weird]);
            id_weird    = unique(id_weird_t);
            
            weird_wb = sort(unique([weird_wb id_weird]))

%             a_dev_theta_L_t2  = [];
%             a_dev_eta_L_t2    = [];
%             a_dev_phi_L_t2    = [];
%             a_dev_theta_R_t2  = [];
%             a_dev_eta_R_t2    = [];
%             a_dev_phi_R_t2    = [];
%
%             for i = 1:nr_wb_t
% 
%                 if isempty(find(id_weird == i,1)) == 1
% 
%                     a_dev_theta_L_t2    = [a_dev_theta_L_t2 a_dev_theta_L_t(:,i)];
%                     a_dev_eta_L_t2      = [a_dev_eta_L_t2 a_dev_eta_L_t(:,i)];
%                     a_dev_phi_L_t2      = [a_dev_phi_L_t2 a_dev_phi_L_t(:,i)];
%                     a_dev_theta_R_t2    = [a_dev_theta_R_t2 a_dev_theta_R_t(:,i)];
%                     a_dev_eta_R_t2      = [a_dev_eta_R_t2 a_dev_eta_R_t(:,i)];
%                     a_dev_phi_R_t2      = [a_dev_phi_R_t2 a_dev_phi_R_t(:,i)];
% 
%                 end
% 
%             end
% 
%             a_dev_theta_L_star   = a_dev_theta_L_t2;
%             a_dev_eta_L_star     = a_dev_eta_L_t2;
%             a_dev_phi_L_star     = a_dev_phi_L_t2;
%             a_dev_theta_R_star   = a_dev_theta_R_t2;
%             a_dev_eta_R_star     = a_dev_eta_R_t2;
%             a_dev_phi_R_star     = a_dev_phi_R_t2;
%             nr_wb_star           = size(a_dev_theta_L_t2,2);
            
        end
        
    end
    

end