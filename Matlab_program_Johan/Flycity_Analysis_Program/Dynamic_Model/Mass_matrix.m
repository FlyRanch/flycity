function [ M_matrix ] = Mass_matrix( body_model, wing_model, body_kin, wing_kin )

    % Determine the mass matrix of the kinetic energy of the fruitfly
    % system and the derivatives of the blocks of the mass matrix.
    
    m_body  = body_model.mass_body;
    m_vw    = wing_model.virtual_mass;
    m_fly   = m_body+2*m_vw;
    
    vb      = body_kin.vb;
    wb      = body_kin.wb;
    
    RL      = wing_kin.RL;
    RR      = wing_kin.RR;
    
    wL      = wing_kin.wL;
    wR      = wing_kin.wR;
    wL_b    = wing_kin.wL_b;
    wR_b    = wing_kin.wR_b;
    
%     R_JL    = zeros(3,1);
%     R_JR    = zeros(3,1);

    cg_b    = body_model.cg_b;
    
    R_JL    = body_model.Joint_left-cg_b;
    R_JR    = body_model.Joint_right-cg_b;
    
    R_cg_L  = wing_model.wing_cg_L;
    R_cg_R  = wing_model.wing_cg_R;
    
%     cg_system = CG_system( m_body, m_vw, R_JL, R_JR, R_cg_L, R_cg_R, RL, RR );
% 
%   cg_b    = cg_system;
%    
%     R_JL    = R_JL-cg_b;
%     R_JR    = R_JR-cg_b;
    
    D_body = [ (cg_b(2)^2+cg_b(3)^2) -cg_b(1)*cg_b(2) -cg_b(1)*cg_b(3); ...
                -cg_b(1)*cg_b(2) (cg_b(1)^2+cg_b(3)^2) -cg_b(2)*cg_b(3); ...
                -cg_b(1)*cg_b(3) -cg_b(2)*cg_b(3) (cg_b(1)^2+cg_b(2)^2)];
      
%     D_L = [ (R_cg_L(2)^2+R_cg_L(3)^2) -R_cg_L(1)*R_cg_L(2) -R_cg_L(1)*R_cg_L(3); ...
%             -R_cg_L(1)*R_cg_L(2) (R_cg_L(1)^2+R_cg_L(3)^2) -R_cg_L(2)*R_cg_L(3); ...
%             -R_cg_L(1)*R_cg_L(3) -R_cg_L(2)*R_cg_L(3) (R_cg_L(1)^2+R_cg_L(2)^2)];
%         
%     D_R = [ (R_cg_R(2)^2+R_cg_R(3)^2) -R_cg_R(1)*R_cg_R(2) -R_cg_R(1)*R_cg_R(3); ...
%             -R_cg_R(1)*R_cg_R(2) (R_cg_R(1)^2+R_cg_R(3)^2) -R_cg_R(2)*R_cg_R(3); ...
%             -R_cg_R(1)*R_cg_R(3) -R_cg_R(2)*R_cg_R(3) (R_cg_R(1)^2+R_cg_R(2)^2)];

    I_body  = body_model.Inertia+m_body*D_body;
    
%     I_vw_L    = wing_model.virtual_Inertia+m_vw*D_L
%     I_vw_R    = wing_model.virtual_Inertia+m_vw*D_R
% 
%     I_body  = body_model.Inertia;

    I_vw      = wing_model.virtual_Inertia;
    
    I_vw_L    = [I_vw(1,1) -I_vw(1,2) I_vw(1,3); ...
                 -I_vw(2,1) I_vw(2,2) I_vw(2,3); ...
                 I_vw(3,1) I_vw(3,2) I_vw(3,3)];
             
    I_vw_R    = I_vw;
    

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
    
    M11     = eye(3)*m_fly;                                                                                                                             % kg
    M12     = m_vw*(-c_JL-RL'*c_cg_L*RL-c_JR-RR'*c_cg_R*RR);                                                                                            % kg*mm
    M13     = m_vw*(-RL'*c_cg_L*RL);                                                                                                                    % kg*mm
    M14     = m_vw*(-RR'*c_cg_R*RR);                                                                                                                    % kg*mm
    M21     = M12';                                                                                                                                     % kg*mm
    M22     = I_body+RL'*I_vw_L*RL+RR'*I_vw_R*RR+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-c_JL-RL'*c_cg_L*RL)+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-c_JR-RR'*c_cg_R*RR);   % kg*mm^2
    M23     = RL'*I_vw_L*RL+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_cg_L*RL);                                                                               % kg*mm^2
    M24     = RR'*I_vw_R*RR+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_cg_R*RR);                                                                               % kg*mm^2
    M31     = M13';                                                                                                                                     % kg*mm
    M32     = M23';                                                                                                                                     % kg*mm^2
    M33     = RL'*I_vw_L*RL+m_vw*(-RL'*c_cg_L*RL)'*(-RL'*c_cg_L*RL);                                                                                    % kg*mm^2
    M34     = zeros(3);                                                                                                                                 % kg*mm^2
    M41     = M14';                                                                                                                                     % kg*mm
    M42     = M24';                                                                                                                                     % kg*mm^2
    M43     = zeros(3);                                                                                                                                 % kg*mm^2
    M44     = RR'*I_vw_R*RR+m_vw*(-RR'*c_cg_R*RR)'*(-RR'*c_cg_R*RR);                                                                                    % kg*mm^2
    
    M       = [ M11 M12 M13 M14; ...
                M21 M22 M23 M24; ...
                M31 M32 M33 M34; ...
                M41 M42 M43 M44];
            
    % Set up the derivatives of the mass matrix over time:

    WL      = [ 0 -wL(3) wL(2); ...
                wL(3) 0 -wL(1); ...
               -wL(2) wL(1) 0];             % rad/s
    
    WR      = [ 0 -wR(3) wR(2); ...
                wR(3) 0 -wR(1); ...
               -wR(2) wR(1) 0];             % rad/s

    WL_R_cg_L = WL*R_cg_L;                  % rad*mm/s
    WR_R_cg_R = WR*R_cg_R;                  % rad*mm/s
    
    c_WL_cg_L = [ 0 -WL_R_cg_L(3) WL_R_cg_L(2); ...
                  WL_R_cg_L(3) 0 -WL_R_cg_L(1); ...
                 -WL_R_cg_L(2) WL_R_cg_L(1) 0]; % rad*mm/s
             
    c_WR_cg_R = [ 0 -WR_R_cg_R(3) WR_R_cg_R(2); ...
                  WR_R_cg_R(3) 0 -WR_R_cg_R(1); ...
                 -WR_R_cg_R(2) WR_R_cg_R(1) 0]; % rad*mm/s

    M12_dot = m_vw*(-RL'*c_WL_cg_L*RL-RR'*c_WR_cg_R*RR);                                        % kg*rad*mm/s
    M13_dot = m_vw*(-RL'*c_WL_cg_L*RL);                                                         % kg*rad*mm/s
    M14_dot = m_vw*(-RR'*c_WR_cg_R*RR);                                                         % kg*rad*mm/s
    M21_dot = M12_dot';                                                                         % kg*rad*mm/s
    M22_dot = RL'*WL*I_vw_L*RL+RR'*WR*I_vw_R*RR+RL'*I_vw_L*WL'*RL+RR'*I_vw_R*WR'*RR+...
              m_vw*(-RL'*c_WL_cg_L*RL)'*(-c_JL-RL'*c_cg_L*RL)+ ...
              m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_WL_cg_L*RL)+ ...
              m_vw*(-RR'*c_WR_cg_R*RR)'*(-c_JR-RR'*c_cg_R*RR)+ ...
              m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_WR_cg_R*RR);                                  % kg*mm^2*rad/s
    M23_dot = RL'*WL*I_vw_L*RL+RL'*I_vw_L*WL'*RL+m_vw*(-RL'*c_WL_cg_L*RL)'* ...
              (-RL'*c_cg_L*RL)+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_WL_cg_L*RL);                 % kg*mm^2*rad/s
    M24_dot = RR'*WR*I_vw_R*RR+RR'*I_vw_R*WR'*RR+m_vw*(-RR'*c_WR_cg_R*RR)'* ...
              (-RR'*c_cg_R*RR)+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_WR_cg_R*RR);                 % kg*mm^2*rad/s

%     M12_dot = m_vw*(-RL'*c_WL_cg_L*RL-RR'*c_WR_cg_R*RR);
%     M13_dot = m_vw*(-RL'*c_WL_cg_L*RL);
%     M14_dot = m_vw*(-RR'*c_WR_cg_R*RR);
%     M21_dot = M12_dot';
%     M22_dot = RL'*WL'*I_vw_L*RL+RR'*WR'*I_vw_R*RR+RL'*I_vw_L*WL*RL+RR'*I_vw_R*WR*RR+...
%               m_vw*(-RL'*c_WL_cg_L*RL)'*(-c_JL-RL'*c_cg_L*RL)+ ...
%               m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_WL_cg_L*RL)+ ...
%               m_vw*(-RR'*c_WR_cg_R*RR)'*(-c_JR-RR'*c_cg_R*RR)+ ...
%               m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_WR_cg_R*RR);
%     M23_dot = RL'*WL'*I_vw_L*RL+RL'*I_vw_L*WL*RL+m_vw*(-RL'*c_WL_cg_L*RL)'* ...
%               (-RL'*c_cg_L*RL)+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_WL_cg_L*RL);
%     M24_dot = RR'*WR'*I_vw_R*RR+RR'*I_vw_R*WR*RR+m_vw*(-RR'*c_WR_cg_R*RR)'* ...
%               (-RR'*c_cg_R*RR)+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_WR_cg_R*RR);


%     M12_dot = m_vw*(-RL'*WL'*c_cg_L*RL-RL'*c_cg_L*WL*RL-RR'*WR'*c_cg_R*RR-RR'*c_cg_R*WR*RR);
%     M13_dot = m_vw*(-RL'*WL'*c_cg_L*RL-RL'*c_cg_L*WL*RL);
%     M14_dot = m_vw*(-RR'*WR'*c_cg_R*RR-RR'*c_cg_R*WR*RR);
%     M21_dot = M12_dot';
%     M22_dot = RL'*WL'*I_vw_L*RL+RR'*WR'*I_vw_R*RR+RL'*I_vw_L*WL*RL+RR'*I_vw_R*WR*RR+...
%               m_vw*(-RL'*WL'*c_cg_L*RL)'*(-c_JL-RL'*c_cg_L*RL)+ ...
%               m_vw*(-RL'*c_cg_L*WL*RL)'*(-c_JL-RL'*c_cg_L*RL)+ ...
%               m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*WL'*c_cg_L*RL)+ ...
%               m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_cg_L*WL*RL)+ ...
%               m_vw*(-RR'*WR'*c_cg_R*RR)'*(-c_JR-RR'*c_cg_R*RR)+ ...
%               m_vw*(-RR'*c_cg_R*WR*RR)'*(-c_JR-RR'*c_cg_R*RR)+ ...
%               m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*WR'*c_cg_R*RR)+ ...
%               m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_cg_R*WR*RR);
%     M23_dot = RL'*WL'*I_vw_L*RL+RL'*I_vw_L*WL*RL+m_vw*(-RL'*WL'*c_cg_L*RL)'* ...
%               (-RL'*c_cg_L*RL)+m_vw*(-RL'*c_cg_L*WL*RL)'* ...
%               (-RL'*c_cg_L*RL)+m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*WL'*c_cg_L*RL)+...
%               +m_vw*(-c_JL-RL'*c_cg_L*RL)'*(-RL'*c_cg_L*WL*RL);
%     M24_dot = RR'*WR'*I_vw_R*RR+RR'*I_vw_R*WR*RR+m_vw*(-RR'*WR'*c_cg_R*RR)'* ...
%               (-RR'*c_cg_R*RR)+m_vw*(-RR'*c_cg_R*WR*RR)'* ...
%               (-RR'*c_cg_R*RR)+m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*WR'*c_cg_R*RR)+...
%               +m_vw*(-c_JR-RR'*c_cg_R*RR)'*(-RR'*c_cg_R*WR*RR);

    % Determine Langrangian derivatives:
    
%     dL_dvb  = 0.5*(vb'*M11+wb'*M21+wL_b'*M31+wR_b'*M41)'+0.5*(M11*vb+M12*wb+M13*wL_b+M14*wR_b);         %kg*mm/s
%     dL_dwb  = 0.5*(vb'*M12+wb'*M22+wL_b'*M32+wR_b'*M42)'+0.5*(M21*vb+M22*wb+M23*wL_b+M24*wR_b);         %kg*mm^2/s

    dL_dvb  = 0.5*(wb'*M21+wL_b'*M31+wR_b'*M41)'+0.5*(M12*wb+M13*wL_b+M14*wR_b);                        %kg*mm/s
    dL_dwb  = 0.5*(vb'*M12+wb'*M22+wL_b'*M32+wR_b'*M42)'+0.5*(M21*vb+M22*wb+M23*wL_b+M24*wR_b);         %kg*mm^2/s

    Vb      = [ 0 -vb(3) vb(2); ...
                vb(3) 0 -vb(1); ...
                -vb(2) vb(1) 0];        % mm/s
    
    Wb      = [ 0 -wb(3) wb(2); ...
                wb(3) 0 -wb(1); ...
                -wb(2) wb(1) 0];        % rad/s
            
    Wb_dL_dvb = Wb*dL_dvb;              % kg*mm/s^2
    Vb_dL_dvb = Vb*dL_dvb;              % kg*mm^2/s^2
    Wb_dL_dwb = Wb*dL_dwb;              % kg*mm^2/s^2
    
    % Determine linear and anguar momentum:
       
    RbwL = R_JL+RL'*R_cg_L;
    RbwR = R_JR+RR'*R_cg_R;
   
    c_bw_L = [0 -RbwL(3) RbwL(2); ...
              RbwL(3) 0 -RbwL(1); ...
              -RbwL(2) RbwL(1) 0];
    c_bw_R = [0 -RbwR(3) RbwR(2); ...
              RbwR(3) 0 -RbwR(1); ...
              -RbwR(2) RbwR(1) 0];
    
%     Mom_mat_b = [ eye(3)*m_fly          -m_vw*(c_bw_L+c_bw_R); ...
%                   m_vw*(c_bw_L+c_bw_R)  I_body+RL'*I_vw_L*RL+RR'*I_vw_R*RR-m_vw*(c_bw_L^2+c_bw_R^2)];

    Mom_mat_b = [ eye(3)*m_fly          -m_vw*(c_bw_L+c_bw_R); ...
                  m_vw*(c_bw_L+c_bw_R)  I_body+RL'*I_vw_L*RL+RR'*I_vw_R*RR-m_vw*(c_bw_L*c_bw_L+c_bw_R*c_bw_R)];
    
    Mom_mat_w = [ m_vw*RL'*c_cg_L*RL    m_vw*RR'*c_cg_R*RR; ...
                  -RL'*I_vw_L*RL+m_vw*c_bw_L*(RL'*c_cg_L*RL)    -RR'*I_vw_R*RR+m_vw*c_bw_R*(RR'*c_cg_R*RR) ];
              
    Momentum = Mom_mat_b*[vb; wb] - Mom_mat_w*[wL_b; wR_b];


    Lin_momentum = Momentum(1:3);
    Ang_momentum = Momentum(4:6);
    
    vb_wb_0 = Mom_mat_b\(Mom_mat_w*[wL_b; wR_b]);
    
    vb_0 = vb_wb_0(1:3);
    wb_0 = vb_wb_0(4:6);
    
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

end

