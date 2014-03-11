function [ alfa_L, beta_L, alfa_R, beta_R ] = wing_alfa_beta( u_wing_L, v_wing_L, w_wing_L, u_wing_R, v_wing_R, w_wing_R )

    % Calculate the angle of attack and angle of sideslip on the wing at
    % different wing locations.
    
    N = length(u_wing_L);
    
    alfa_L = zeros(1,N);
    
    beta_L = zeros(1,N);
    
    alfa_R = zeros(1,N);
    
    beta_R = zeros(1,N);
    
    for i = 1:N
        
        alfa_L(i) = real(atan2(w_wing_L(i),u_wing_L(i)));
        beta_L(i) = real(atan2(v_wing_L(i),u_wing_L(i)));
        alfa_R(i) = real(atan2(w_wing_R(i),u_wing_R(i)));
        beta_R(i) = real(atan2(v_wing_R(i),u_wing_R(i)));
        
    end
    
end

