function [ Fg_b, Fg_strkpln ] = Gravity_instantaneous( body_model, kine )


    % Compute the direction and magnitude of the gravity vector in the body
    % reference frame:
    
    m_fly = body_model.mass_fly;
    
    Rb      = kine.Rb;
    R_strk  = kine.R_strk;
    g       = -body_model.g; %[m/s^2]
    
    Fg_world = g*[0; 0; m_fly]; % kg * [m/s^2] = N
    
    Fg_b = Rb*Fg_world;
    Fg_strkpln = R_strk*Fg_b;

end

