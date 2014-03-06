function [ gravity_f ] = Gravity_forces( q_b , mass_fly, n )


    % Compute direction gravity vector for the orientation of the
    % strokeplane:
    
    gravity_f = {};
    
    beta = (55/180)*pi;
    
    R_beta = [cos(-beta) 0 -sin(-beta); ...
              0 1 0; ...
              sin(-beta) 0 cos(-beta)]; 
          
    Fg_b = zeros(3,n);
    
    Fg_st = zeros(3,n);
    
    R_b = zeros(3,3,n);
    
    for i = 1:n
        
        R_b(:,:,i) = quat2matNEW([q_b(1,i); q_b(2,i); q_b(3,i); q_b(4,i)]);
        
        Fg_b(:,i) = R_b(:,:,i)'*[0; 0; -9810*mass_fly]; % [kg*mm/s^2]
        
        Fg_st(:,i) = R_beta*Fg_b(:,i);
        
    end
    
    gravity_f.R_b = R_b;
    gravity_f.Fg_b = Fg_b;
    gravity_f.Fg_st = Fg_st;
    


end

