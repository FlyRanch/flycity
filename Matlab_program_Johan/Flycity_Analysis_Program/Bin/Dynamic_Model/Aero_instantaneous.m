function [ FM_b, FM_strkpln, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aero_instantaneous( kine, body_model, wing_model, rot_lift_on )

    % Compute the aerodynamic forces on the wings by means of a
    % quasi-steady model:
    
    vb              = kine.vb;
    wb              = kine.wb;
    wL              = kine.wL;              % [ rad/s ]
    wR              = kine.wR;              % [ rad/s ]
    w_dot_L         = kine.w_dot_L;
    w_dot_R         = kine.w_dot_R;
    RL              = kine.RL;
    RR              = kine.RR;
    R_strk          = kine.R_strk;
    
    Joint_left      = body_model.Joint_left;    % [ mm ]
    Joint_right     = body_model.Joint_right;   % [ mm ]
    cg_b            = body_model.cg_b;          % [ mm ]
    
    y_sect_L        = wing_model.y_sect_L;      % [ mm ]
    chords_L        = wing_model.chords_L;      % [ mm ]
    x_LE_L          = wing_model.x_LE_L;        % [ mm ]
    y_sect_R        = wing_model.y_sect_R;      % [ mm ]
    chords_R        = wing_model.chords_R;      % [ mm ]
    x_LE_R          = wing_model.x_LE_R;        % [ mm ]
    rho             = wing_model.rho;           % [ kg/mm^3 ]
    
    delta_R         = abs(y_sect_L(2,2)-y_sect_L(2,1));     % [ mm ]
    
    C_rot           = 1.55;
    
    nr_sect = length(y_sect_L);
    
    FM_strkpln  = nan(6,1);

    % Construct U, U_dot, alfa and alfa_dot

    UL          = nan(3,nr_sect);                   % [ mm/s ]
    UR          = nan(3,nr_sect);                   % [ mm/s ]
    U_dot_L     = nan(3,nr_sect);                   % [mm/s^2]
    U_dot_R     = nan(3,nr_sect);                   % [mm/s^2]
    alfa_L      = nan(nr_sect,1);                   % [ rad ]
    alfa_R      = nan(nr_sect,1);                   % [ rad ]
    alfa_dot_L  = nan(nr_sect,1);                   % [ rad/s ]
    alfa_dot_R  = nan(nr_sect,1);                   % [ rad/s ]
    
    for j = 1:nr_sect
                
        % Compute the local left and right wing velocity at the 
        % current section:

        UL(:,j) = cross(wL,y_sect_L(:,j))+cross(RL*wb,RL*Joint_left+y_sect_L(:,j))+RL*vb;   % [ mm/s ]
        UR(:,j) = cross(wR,y_sect_R(:,j))+cross(RR*wb,RR*Joint_right+y_sect_R(:,j))+RR*vb;  % [ mm/s ]
        
%         WL      = [ 0 -wL(3) wL(2); ...
%                     wL(3) 0 -wL(1); ...
%                    -wL(2) wL(1) 0];             % rad/s
% 
%         WR      = [ 0 -wR(3) wR(2); ...
%                     wR(3) 0 -wR(1); ...
%                    -wR(2) wR(1) 0];             % rad/s
%         
%         U_dot_L(:,j) = cross(w_dot_L,y_sect_L(:,j))+cross(WL*RL*wb,RL*Joint_left+y_sect_L(:,j))+...
%                        cross(RL*wb,WL*RL*Joint_left+y_sect_L(:,j))+WL*RL*vb;
%         
%         U_dot_R(:,j) = cross(w_dot_R,y_sect_R(:,j))+cross(WR*RR*wb,RR*Joint_right+y_sect_R(:,j))+...
%                        cross(RR*wb,WR*RR*Joint_right+y_sect_R(:,j))+WR*RR*vb;

        alfa_L(j) = real(atan2(UL(3,j),UL(1,j)));     % [ rad ]
        alfa_R(j) = real(atan2(UR(3,j),UR(1,j)));     % [ rad ]
%         alfa_dot_L(j) = -(UL(3,j)/(UL(1,j)^2+UL(3,j)^2))*U_dot_L(1,j)+(UL(1,j)/(UL(1,j)^2+UL(3,j)^2))*U_dot_L(3,j); % [ rad/s ]
%         alfa_dot_R(j) = -(UR(3,j)/(UR(1,j)^2+UR(3,j)^2))*U_dot_R(1,j)+(UR(1,j)/(UR(1,j)^2+UR(3,j)^2))*U_dot_R(3,j); % [ rad/s ]
                
        alfa_dot_L(j) = -wL(2);
        alfa_dot_R(j) = -wR(2);

    end
    
    F_trans_L = zeros(nr_sect,3); % [ N ]
    F_trans_R = zeros(nr_sect,3); % [ N ]
    M_trans_L = zeros(nr_sect,3); % [ N*mm ]
    M_trans_R = zeros(nr_sect,3); % [ N*mm ]
    F_rot_L   = zeros(nr_sect,3); % [ N ]
    F_rot_R   = zeros(nr_sect,3); % [ N ]
    M_rot_L   = zeros(nr_sect,3); % [ N*mm ]
    M_rot_R   = zeros(nr_sect,3); % [ N*mm ]

    
    for j = 1:nr_sect
        
        % Compute location center of pressure:
            
        x_cp_L = x_LE_L(j)-chords_L(j)*(0.82*abs(alfa_L(j))/pi+0.05);
        x_cp_R = x_LE_R(j)-chords_R(j)*(0.82*abs(alfa_R(j))/pi+0.05);
    
        if alfa_L(j) >= 0
            
            alfa_left        = abs(radtodeg(alfa_L(j)));
            Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
            Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
            L_trans_L        = 1e-3*0.5*rho*chords_L(j)*(UL(1,j)^2+UL(3,j)^2)*Cl_left*delta_R;        % [ N ]
            D_trans_L        = 1e-3*0.5*rho*chords_L(j)*(UL(1,j)^2+UL(3,j)^2)*Cd_left*delta_R;        % [ N ]
            F_trans_L(j,:)   = [(L_trans_L*sind(alfa_left)-D_trans_L*cosd(alfa_left)) ...
                                0 ...
                                -(L_trans_L*cosd(alfa_left)+D_trans_L*sind(alfa_left))];                  % [ N ]            
            M_trans_L(j,:)   = (cross([x_cp_L; y_sect_L(2,j); y_sect_L(3,j)],F_trans_L(j,:)'))';
            
        else
            
            alfa_left        = abs(radtodeg(alfa_L(j)));
            Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
            Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
            L_trans_L        = 1e-3*0.5*rho*chords_L(j)*(UL(1,j)^2+UL(3,j)^2)*Cl_left*delta_R;        % [ N ]
            D_trans_L        = 1e-3*0.5*rho*chords_L(j)*(UL(1,j)^2+UL(3,j)^2)*Cd_left*delta_R;        % [ N ]
            F_trans_L(j,:)   = [(L_trans_L*sind(alfa_left)-D_trans_L*cosd(alfa_left)) ...
                                0 ...
                                (L_trans_L*cosd(alfa_left)+D_trans_L*sind(alfa_left))];                  % [ N ]          
            M_trans_L(j,:)   = (cross([x_cp_L; y_sect_L(2,j); y_sect_L(3,j)],F_trans_L(j,:)'))';
            
        end
        
        if rot_lift_on == 1
                            
            F_rot_L(j,:) = 1e-3*[ 0 0 C_rot*rho*alfa_dot_L(j)*chords_L(j)^2*sqrt(UL(1,j)^2+UL(3,j)^2)*delta_R];     % [ N ]
            M_rot_L(j,:) = (cross([x_cp_L; y_sect_L(2,j); y_sect_L(3,j)],F_rot_L(j,:)'))';                         % [ N*mm ]            
                            
        end   
        
        if alfa_R(j) >= 0
            
            alfa_right       = abs(radtodeg(alfa_R(j)));
            Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
            Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
            L_trans_R        = 1e-3*0.5*rho*chords_R(j)*(UR(1,j)^2+UR(3,j)^2).*Cl_right*delta_R;       % [ N ]
            D_trans_R        = 1e-3*0.5*rho*chords_R(j)*(UR(1,j)^2+UR(3,j)^2).*Cd_right*delta_R;       % [ N ]
            F_trans_R(j,:)   = [(L_trans_R*sind(alfa_right)-D_trans_R*cosd(alfa_right)) ...
                                0 ...
                                -(L_trans_R*cosd(alfa_right)+D_trans_R*sind(alfa_right))];             % [ N ]
            M_trans_R(j,:)   = (cross([x_cp_R; y_sect_R(2,j); y_sect_R(3,j)],F_trans_R(j,:)'))';
                
        else
            
            alfa_right       = abs(radtodeg(alfa_R(j)));
            Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
            Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
            L_trans_R        = 1e-3*0.5*rho*chords_R(j)*(UR(1,j)^2+UR(3,j)^2).*Cl_right*delta_R;       % [ N ]
            D_trans_R        = 1e-3*0.5*rho*chords_R(j)*(UR(1,j)^2+UR(3,j)^2).*Cd_right*delta_R;       % [ N ]
            F_trans_R(j,:)   = [(L_trans_R*sind(alfa_right)-D_trans_R*cosd(alfa_right)) ...
                                0 ...
                                (L_trans_R*cosd(alfa_right)+D_trans_R*sind(alfa_right))];             % [ N ]
            M_trans_R(j,:)   = (cross([x_cp_R; y_sect_R(2,j); y_sect_R(3,j)],F_trans_R(j,:)'))';

        end
        
        if rot_lift_on == 1
                            
            F_rot_R(j,:) = 1e-3*[ 0 0 C_rot*rho*alfa_dot_R(j)*chords_R(j)^2*sqrt(UR(1,j)^2+UR(3,j)^2)*delta_R];     % [ N ]
            M_rot_R(j,:) = (cross([x_cp_R; y_sect_R(2,j); y_sect_R(3,j)],F_rot_R(j,:)'))';                         % [ N*mm ]            
                            
        end
        
    end
            
    FM_L = [sum(F_trans_L)'+sum(F_rot_L)'; sum(M_trans_L)'+sum(M_rot_L)'];
    FM_R = [sum(F_trans_R)'+sum(F_rot_R)'; sum(M_trans_R)'+sum(M_rot_R)'];
    
    FM_b = nan(6,1);
    
    FM_b(1:3)       = RL'*FM_L(1:3)+RR'*FM_R(1:3);
    FM_b(4:6)       = RL'*FM_L(4:6)+RR'*FM_R(4:6)+cross((Joint_left-cg_b),RL'*FM_L(1:3))+cross((Joint_right-cg_b),RR'*FM_R(1:3));

%     FM_b(1:3)       = RL'*FM_L(1:3)+RR'*FM_R(1:3);
%     FM_b(4:6)       = RL'*FM_L(4:6)+RR'*FM_R(4:6);

    FM_strkpln(1:3) = R_strk*FM_b(1:3);
    FM_strkpln(4:6) = R_strk*FM_b(4:6);


end

