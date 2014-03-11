function [ M_matrix ] = Mass_matrix( body_model, wing_model, body_kin, wing_kin )

    % Determine the mass matrix of the kinetic energy of the fruitfly
    % system and the derivatives of the blocks of the mass matrix.
    
    % Obtain body and wing parameters and body and wing kinematics:
    
    m_body  = body_model.mass_body;
    m_vw    = wing_model.virtual_mass;
    m_fly   = m_body+2*m_vw;
    
    vb      = body_kin.vb;
    wb      = body_kin.wb;
    
    RL      = wing_kin.RL;
    RR      = wing_kin.RR;
    
    wL_b    = wing_kin.wL_b;
    wR_b    = wing_kin.wR_b;
    
    cg_b    = body_model.cg_b;
    
    R_JL    = body_model.Joint_left-cg_b;
    R_JR    = body_model.Joint_right-cg_b;
    R_cg_L  = wing_model.wing_cg_L;
    R_cg_R  = wing_model.wing_cg_R;

        
    I_body  = body_model.Inertia-m_body*(cg_b*cg_b');
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
           
    % Set up the mass matrix:
    
    M11     = eye(3)*m_fly;
    M12     = m_vw*(-c_JL-RL'*c_cg_L*RL-c_JR-RR'*c_cg_R*RR);
    M13     = m_vw*(-RL'*c_cg_L*RL);
    M14     = m_vw*(-RR'*c_cg_R*RR);
    M21     = M12';
    M22     = I_body+RL'*I_vw*RL+RR'*I_vw*RR+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-c_JL-RL'*c_cg_L*RL)+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-c_JR-RR'*c_cg_R*RR);
    M23     = RL'*I_vw*RL+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_cg_L*RL);
    M24     = RR'*I_vw*RR+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_cg_R*RR);
    M31     = M13';
    M32     = M23';
    M33     = RL'*I_vw*RL+m_vw*(-RL'*c_cg_L*RL)'*(-RL'*c_cg_L*RL);
    M34     = zeros(3);
    M41     = M14';
    M42     = M24';
    M43     = zeros(3);
    M44     = RR'*I_vw*RR+m_vw*(-RR'*c_cg_R*RR)'*(-RR'*c_cg_R*RR);
    
    M       = [ M11 M12 M13 M14; ...
                M21 M22 M23 M24; ...
                M31 M32 M33 M34; ...
                M41 M42 M43 M44];
            
    % Set up the derivatives of the mass matrix over time:
    
    WL      = [ 0 -wL_b(3) wL_b(2); ...
                wL_b(3) 0 -wL_b(1); ...
               -wL_b(2) wL_b(1) 0];
    
    WR      = [ 0 -wR_b(3) wR_b(2); ...
                wR_b(3) 0 -wR_b(1); ...
               -wR_b(2) wR_b(1) 0];
    
    M12_dot = m_vw*(-RL'*WL'*c_cg_L*RL-RL'*c_cg_L*WL*RL-RR'*WR'*c_cg_R*RR-RR'*c_cg_R*WR*RR);
    M13_dot = m_vw*(-RL'*WL'*c_cg_L*RL-RL'*c_cg_L*WL*RL);
    M14_dot = m_vw*(-RR'*WR'*c_cg_R*RR-RR'*c_cg_R*WR*RR);
    M21_dot = M12_dot';
    M22_dot = RL'*WL'*I_vw*RL+RL'*I_vw*WL*RL+RR'*WR'*I_vw*RR+RR'*I_vw*WR*RR+ ...
              m_vw*(-RL'*WL'*c_JL*RL-RL'*c_JL*WL*RL)'*(-c_JL-RL'*c_cg_L*RL)+ ...
              m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*WL'*c_JL*RL-RL'*c_JL*WL*RL)+ ...
              m_vw*(-RR'*WR'*c_JR*RR-RR'*c_JR*WR*RR)'*(-c_JR-RR'*c_cg_R*RR)+ ...
              m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*WR'*c_JR*RR-RR'*c_JR*WR*RR);
    M23_dot = RL'*WL'*I_vw*RL+RL'*I_vw*WL*RL+m_vw*(-RL'*WL'*c_cg_L*RL-RL'*c_cg_L*WL*RL)'* ...
              (-RL'*c_cg_L*RL)+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*WL'*c_cg_L*RL-RL'*c_cg_L*WL*RL);
    M24_dot = RR'*WR'*I_vw*RR+RR'*I_vw*WR*RR+m_vw*(-RR'*WR'*c_cg_R*RR-RR'*c_cg_R*WR*RR)'* ...
              (-RR'*c_cg_R*RR)+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*WR'*c_cg_R*RR-RR'*c_cg_R*WR*RR);
          
    % Determine Langrangian derivatives:
    
    dL_dvb  = 0.5*(vb'*M11+wb'*M21+wL_b'*M31+wR_b'*M41)'+0.5*(M11*vb+M12*wb+M13*wL_b+M14*wR_b);
    dL_dwb  = 0.5*(vb'*M12+wb'*M22+wL_b'*M32+wR_b'*M42)'+0.5*(M21*vb+M22*wb+M23*wL_b+M24*wR_b);
    
    Vb      = [ 0 -vb(3) vb(2); ...
                vb(3) 0 -vb(1); ...
                -vb(2) vb(1) 0];
    
    Wb      = [ 0 -wb(3) wb(2); ...
                wb(3) 0 -wb(1); ...
                -wb(2) wb(1) 0];
            
    Wb_dL_dvb = Wb*dL_dvb;
    Vb_dL_dvb = Vb*dL_dvb;
    Wb_dL_dwb = Wb*dL_dwb;
    
    % Determine linear and anguar momentum:
    
    RbwL = R_JL+RL'*R_cg_L;
    RbwR = R_JR+RR'*R_cg_R;
    
    c_bw_L = [0 -RbwL(3) RbwL(2); ...
              RbwL(3) 0 -RbwL(1); ...
              -RbwL(2) RbwL(1) 0];
    c_bw_R = [0 -RbwR(3) RbwR(2); ...
              RbwR(3) 0 -RbwR(1); ...
              -RbwR(2) RbwR(1) 0];
    
    Mom_mat_b = [ eye(3)*m_fly          -m_vw*(c_bw_L+c_bw_R); ...
                  m_vw*(c_bw_L+c_bw_R)  I_body+RL'*I_vw*RL+RR'*I_vw*RR-m_vw*(c_bw_L^2+c_bw_R^2)];
    
    Mom_mat_w = [ m_vw*RL'*c_cg_L*RL    m_vw*RR'*c_cg_R*RR; ...
                  -RL'*I_vw*RL+m_vw*c_bw_L*RL'*c_cg_L*RL    -RR'*I_vw*RR+m_vw*c_bw_R*RR'*c_cg_R*RR ];
              
    Momentum = Mom_mat_b*[vb; wb] - Mom_mat_w*[wL_b; wR_b];
    
    Lin_momentum = Momentum(1:3);
    Ang_momentum = Momentum(4:6);
    
    vb_wb_0 = Mom_mat_b\(Mom_mat_w*[wL_b; wR_b]);
    
    vb_0 = vb_wb_0(1:3);
    wb_0 = vb_wb_0(4:6);
    
    RcgLwL  = RL'*c_cg_L*RL*wL_b;
    RcgRwR  = RR'*c_cg_R*RR*wR_b;

    % Return the data:
    
    M_matrix = {};
    
    M_matrix.M                  = M;
    
    M_matrix.M11                = M11;
    M_matrix.M12                = M12;
    M_matrix.M13                = M13;
    M_matrix.M14                = M14;
    M_matrix.M21                = M21;
    M_matrix.M22                = M22;
    M_matrix.M23                = M23;
    M_matrix.M24                = M24;
    M_matrix.M31                = M31;
    M_matrix.M32                = M32;
    M_matrix.M33                = M33;
    M_matrix.M34                = M34;
    M_matrix.M41                = M41;
    M_matrix.M42                = M42;
    M_matrix.M43                = M43;
    M_matrix.M44                = M44;
    
    M_matrix.M12_dot            = M12_dot;
    M_matrix.M13_dot            = M13_dot;
    M_matrix.M14_dot            = M14_dot;
    M_matrix.M21_dot            = M21_dot;
    M_matrix.M22_dot            = M22_dot;
    M_matrix.M23_dot            = M23_dot;
    M_matrix.M24_dot            = M24_dot;
    
    M_matrix.Wb_dL_dvb          = Wb_dL_dvb;
    M_matrix.Vb_dL_dvb          = Vb_dL_dvb;
    M_matrix.Wb_dL_dwb          = Wb_dL_dwb;
    
    M_matrix.Lin_momentum       = Lin_momentum;
    M_matrix.Ang_momentum       = Ang_momentum;
    
    M_matrix.vb_0               = vb_0;
    M_matrix.wb_0               = wb_0;
    
    M_matrix.RcgLwL             = RcgLwL;
    M_matrix.RcgRwR             = RcgRwR;
    
end

