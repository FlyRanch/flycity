function [ q_inv ] = quat_inv( q )


    % Conjugate:
    
    q_conj = [-q(1); -q(2); -q(3); q(4)];
    
    % Inverse:
    
    q_inv = q_conj./norm(q)^2;
    
    q_inv = q_inv./norm(q_inv);


end

