for i = 1:seqs
    
    A = trims(i,2); % min
    B = trims(i,3); % max
    C = trims(i,4); % min 2
    D = trims(i,5); % max 2
    
    for j = 1:length(pathDB.phi_R(:,i))
        
        if j > A && j < B
            
            pathDB.phi_R(j,i) = nan;
            pathDB.phi_L(j,i) = nan;
            pathDB.eta_R(j,i) = nan;
            pathDB.eta_L(j,i) = nan;
            pathDB.theta_R(j,i) = nan;
            pathDB.theta_L(j,i) = nan;
            
        elseif j > C && j < D
            
            pathDB.phi_R(j,i) = nan;
            pathDB.phi_L(j,i) = nan;
            pathDB.eta_R(j,i) = nan;
            pathDB.eta_L(j,i) = nan;
            pathDB.theta_R(j,i) = nan;
            pathDB.theta_L(j,i) = nan;
            
        end
    end
end

clear A B C D i j
    