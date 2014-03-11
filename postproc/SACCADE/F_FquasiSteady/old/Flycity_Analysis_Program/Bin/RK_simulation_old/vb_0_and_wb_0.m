function [ vb_0, wb_0 ] = vb_0_and_wb_0( body_model, wing_model, wing_kin )


    % Obtain body and wing parameters and body and wing kinematics:
    
        
    m_body  = body_model.mass_body;
    m_vw    = wing_model.virtual_mass;
    
    wL      = wing_kin.wL;
    wR      = wing_kin.wR;
    
    
    R_JL    = body_model.Joint_left;
    R_JR    = body_model.Joint_right;
    R_cg_L  = wing_model.wing_cg_L;
    R_cg_R  = wing_model.wing_cg_R;
        
    I_body  = body_model.Inertia;
    I_vw    = wing_model.virtual_Inertia;
    
    % Define cross matrices:
    
    c_JL    = [ 0 -R_JL(3) R_JL(2); ...
                R_JL(3) 0 -R_JL(1); ...
               -R_JL(2) R_JL(1) 0];
    
    c_JR    = [ 0 -R_JR(3) R_JR(2); ...
                R_JR(3) 0 -R_JR(1); ...
               -R_JR(2) R_JR(1) 0];
    
    c_cg_L  = [ 0 -R_cg_L(3) R_cg_L(2); ...
                R_cg_L(3) 0 -R_cg_L(1); ...
               -R_cg_L(2) R_cg_L(1) 0];
    
    c_cg_R  = [ 0 -R_cg_R(3) R_cg_R(2); ...
                R_cg_R(3) 0 -R_cg_R(1); ...
               -R_cg_R(2) R_cg_R(1) 0];
    
    % Compute the required body velocity and body angular velocity for zero
    % angular momentum:
    
    N = length(wL(1,:));
    
    vb_0 = nan(3,N);
    wb_0 = nan(3,N);
    
    for i = 1:N
    
    A = [ eye(3)*(m_body+2*m_vw)             (c_JL+c_JR)*m_vw; ...
          (c_JL+c_JR)*m_vw      I_body+2*I_vw-m_vw*(c_JL*c_JL+c_JR*c_JR)];
      
    b = [ m_vw*(c_cg_L*wL(:,i)+c_cg_R*wR(:,i)); ...
          -I_vw*wL(:,i)-I_vw*wR(:,i)+m_vw*(c_JL*c_cg_L*wL(:,i)+c_JR*c_cg_R*wR(:,i))];
      
    c = A\b;
    
    vb_0(:,i) = c(1:3);
    wb_0(:,i) = c(4:6);
    
    end

end

