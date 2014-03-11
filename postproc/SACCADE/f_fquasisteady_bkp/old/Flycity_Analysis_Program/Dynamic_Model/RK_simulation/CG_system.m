function [ cg_system ] = CG_system( mb, m_vw, R_JL, R_JR, R_cg_L, R_cg_R, RL, RR )

    % Compute the location of the center of gravity of the wing-body
    % system. The output vector is in the body reference frame w.r.t. the
    % 0-location of the body.
    
    R_wing_L = R_JL + RL'*R_cg_L;
    R_wing_R = R_JR + RR'*R_cg_R;
    
    % Non-dimensionalize axis of rotation:
             
    % Obtain center of gravity of the system:
    
    cg_system =(1/mb)*(-m_vw*(R_wing_L+R_wing_R));
    
end

