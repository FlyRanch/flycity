function [ FM_b, FM_strkpln, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aero_instantaneous( kine, body_model, wing_model, wingbeat, rot_lift_on )

    % Compute the aerodynamic forces on the wings by means of a
    % quasi-steady model:
    
    vb              = kine.vb;
    wb              = kine.wb;
    wL              = kine.wL;              % [ rad/s ]
    wR              = kine.wR;              % [ rad/s ]
    RL              = kine.RL;
    RR              = kine.RR;
    R_strk          = kine.R_strk;
    alfa_L_old      = kine.alfa_L_old;
    alfa_R_old      = kine.alfa_R_old;
    dt              = kine.dt;
    id              = kine.id;
    
    Joint_left      = body_model.Joint_left;    % [ mm ]
    Joint_right     = body_model.Joint_right;   % [ mm ]
    cg_b            = body_model.cg_b;          % [ mm ]
    
    y_sect_L        = wing_model.y_sect_L;      % [ mm ]
    chords_L        = wing_model.chords_L;      % [ mm ]
    y_sect_R        = wing_model.y_sect_R;      % [ mm ]
    chords_R        = wing_model.chords_R;      % [ mm ]
    rho             = wing_model.rho;           % [ kg/mm^3 ]
    
    wb_loc          = wingbeat.wb_loc;
    down_loc        = wingbeat.down_loc;
    up_loc          = wingbeat.up_loc;
    
    delta_R         = abs(y_sect_L(2,2)-y_sect_L(2,1));     % [ mm ]
    
    C_rot           = 1.55;
    
    nr_sect = length(y_sect_L);
    
    FM_strkpln  = nan(6,1);
    FM_L        = nan(6,1);
    FM_R        = nan(6,1);
    
    UL          = nan(3,nr_sect);                 % [ mm/s ]
    UR          = nan(3,nr_sect);                 % [ mm/s ]
    alfa_L      = nan(nr_sect,1);                   % [ rad ]
    alfa_R      = nan(nr_sect,1);                   % [ rad ]
    alfa_dot_L  = nan(nr_sect,1);                   % [ rad/s ]
    alfa_dot_R  = nan(nr_sect,1);                   % [ rad/s ]

            
    for j = 1:nr_sect
                
        % Compute the local left and right wing velocity at the 
        % current section:
        

        UL(:,j) = cross(wL,y_sect_L(:,j))+cross(RL*wb,RL*Joint_left+y_sect_L(:,j))+RL*vb;   % [ mm/s ]
        UR(:,j) = cross(wR,y_sect_R(:,j))+cross(RR*wb,RR*Joint_right+y_sect_R(:,j))+RR*vb;  % [ mm/s ]

        alfa_L(j) = real(atan2(UL(3,j),UL(1,j)));     % [ rad ]
        alfa_R(j) = real(atan2(UR(3,j),UR(1,j)));     % [ rad ]
                
    end

    alfa_dot_L_t = alfa_L-alfa_L_old;
    alfa_dot_R_t = alfa_R-alfa_R_old;
    
    for k = 1:nr_sect
        
        if isnan(alfa_L_old(k))==1
            
            alfa_dot_L(k) = 0;
            
        else
            
            if abs(alfa_dot_L_t(k)) >= (pi/2)
            
                alfa_dot_L(k) = 0;
            
            else
                
                alfa_dot_L(k) = (alfa_L(k)-alfa_L_old(k))/dt;
            
            end
            
        end
    
        if isnan(alfa_R_old(k))==1
            
            alfa_dot_R(k) = 0;
            
        else
            
            if abs(alfa_dot_R_t(k)) >= (pi/2)
            
                alfa_dot_R(k) = 0;
                
            else
                
                alfa_dot_R(k) = (alfa_R(k)-alfa_R_old(k))/dt;
            
            end
            
        end
        
    end
    
    
%     L_trans_L  = zeros(nr_sect,1);  % [ N ]
%     L_trans_R  = zeros(nr_sect,1);  % [ N ]
%     D_trans_L  = zeros(nr_sect,1);  % [ N ]
%     D_trans_R  = zeros(nr_sect,1);  % [ N ]
    F_trans_L = zeros(nr_sect,3); % [ N ]
    F_trans_R = zeros(nr_sect,3); % [ N ]
    M_trans_L = zeros(nr_sect,3); % [ N*mm ]
    M_trans_R = zeros(nr_sect,3); % [ N*mm ]
    F_rot_L   = zeros(nr_sect,3); % [ N ]
    F_rot_R   = zeros(nr_sect,3); % [ N ]
    M_rot_L   = zeros(nr_sect,3); % [ N*mm ]
    M_rot_R   = zeros(nr_sect,3); % [ N*mm ]
    
    if id >= down_loc(1) && id <= down_loc(2)
                
        alfa_left        = radtodeg(alfa_L);                                                       % [ deg ]
        alfa_right       = radtodeg(alfa_R);                                                       % [ deg ]
        Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
        Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
        Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
        Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
%         L_trans_L        = 1e-3*0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cl_left*delta_R;        % [ N ]
%         L_trans_R        = 1e-3*0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cl_right*delta_R;       % [ N ]
%         D_trans_L        = 1e-3*0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cd_left*delta_R;        % [ N ]
%         D_trans_R        = 1e-3*0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cd_right*delta_R;       % [ N ]
        L_trans_L        = 0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cl_left*delta_R;        % [ N ]
        L_trans_R        = 0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cl_right*delta_R;       % [ N ]
        D_trans_L        = 0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cd_left*delta_R;        % [ N ]
        D_trans_R        = 0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cd_right*delta_R;       % [ N ]
        F_trans_L        = [(L_trans_L.*sind(alfa_left)-D_trans_L.*cosd(alfa_left)) ...
                            zeros(nr_sect,1) ...
                            -(L_trans_L.*cosd(alfa_left)+D_trans_L.*sind(alfa_left))];        % [ N ]
        F_trans_R        = [(L_trans_R.*sind(alfa_right)-D_trans_R.*cosd(alfa_right)) ...
                            zeros(nr_sect,1) ...
                            -(L_trans_R.*cosd(alfa_right)+D_trans_R.*sind(alfa_right))];      % [ N ]

        for j = 1:nr_sect

            M_trans_L(j,:) = (cross(y_sect_L(:,j),F_trans_L(j,:)'))';                     % [ N*mm ]
            M_trans_R(j,:) = (cross(y_sect_R(:,j),F_trans_R(j,:)'))';                     % [ N*mm ]

        end
                
        if rot_lift_on == 1
                    
%             F_rot_L = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_L.*chords_L.^2.*sqrt((UL(1,:)').^2+(UL(3,:)').^2)*delta_R];     % [ N ]
%             F_rot_R = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_R.*chords_R.^2.*sqrt((UR(1,:)').^2+(UR(3,:)').^2)*delta_R];     % [ N ]

            F_rot_L = [ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_L.*chords_L.^2.*sqrt((UL(1,:)').^2+(UL(3,:)').^2)*delta_R];     % [ N ]
            F_rot_R = [ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_R.*chords_R.^2.*sqrt((UR(1,:)').^2+(UR(3,:)').^2)*delta_R];     % [ N ]

            for j = 1:nr_sect

                M_rot_L(j,:) = (cross(y_sect_L(:,j),F_rot_L(j,:)'))';                         % [ N*mm ]
                M_rot_R(j,:) = (cross(y_sect_R(:,j),F_rot_R(j,:)'))';                         % [ N*mm ]

            end
                    
        end
                
    elseif id >= up_loc(1) && id <= up_loc(2)
                
            alfa_left        = -radtodeg(alfa_L);                                                      % [ N ] 
            alfa_right       = -radtodeg(alfa_R);                                                      % [ N ] 
            Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
            Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
            Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
            Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
%             L_trans_L        = 1e-3*0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cl_left*delta_R;        % [ N ] 
%             L_trans_R        = 1e-3*0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cl_right*delta_R;       % [ N ] 
%             D_trans_L        = 1e-3*0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cd_left*delta_R;        % [ N ] 
%             D_trans_R        = 1e-3*0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cd_right*delta_R;       % [ N ] 
            L_trans_L        = 0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cl_left*delta_R;        % [ N ] 
            L_trans_R        = 0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cl_right*delta_R;       % [ N ] 
            D_trans_L        = 0.5*rho*chords_L.*((UL(1,:)').^2+(UL(3,:)').^2).*Cd_left*delta_R;        % [ N ] 
            D_trans_R        = 0.5*rho*chords_R.*((UR(1,:)').^2+(UR(3,:)').^2).*Cd_right*delta_R;       % [ N ] 
            F_trans_L        = [L_trans_L.*sind(alfa_left)-D_trans_L.*cosd(alfa_left) ...
                                zeros(nr_sect,1) ...
                                L_trans_L.*cosd(alfa_left)+D_trans_L.*sind(alfa_left)];          % [ N ] 
            F_trans_R        = [L_trans_R.*sind(alfa_right)-D_trans_R.*cosd(alfa_right) ...
                                zeros(nr_sect,1) ...
                                L_trans_R.*cosd(alfa_right)+D_trans_R.*sind(alfa_right)];        % [ N ] 
                                 
            for j = 1:nr_sect

                M_trans_L(j,:) = (cross(y_sect_L(:,j),F_trans_L(j,:)'))';                     % [ N*mm ]
                M_trans_R(j,:) = (cross(y_sect_R(:,j),F_trans_R(j,:)'))';                     % [ N*mm ]

            end

            if rot_lift_on == 1
                    
%                 F_rot_L = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_L.*chords_L.^2.*sqrt((UL(1,:)').^2+(UL(3,:)').^2)*delta_R];    % [ N ] 
%                 F_rot_R = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_R.*chords_R.^2.*sqrt((UR(1,:)').^2+(UR(3,:)').^2)*delta_R];    % [ N ]

                F_rot_L = [ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_L.*chords_L.^2.*sqrt((UL(1,:)').^2+(UL(3,:)').^2)*delta_R];    % [ N ] 
                F_rot_R = [ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_R.*chords_R.^2.*sqrt((UR(1,:)').^2+(UR(3,:)').^2)*delta_R];    % [ N ]

                for j = 1:nr_sect

                    M_rot_L(j,:) = (cross(y_sect_L(:,j),F_rot_L(j,:)'))';                         % [ N*mm ]
                    M_rot_R(j,:) = (cross(y_sect_R(:,j),F_rot_R(j,:)'))';                         % [ N*mm ]

                end
                    
            end
                
    end
            
    FM_L = 1e-3*[sum(F_trans_L)'+sum(F_rot_L)'; sum(M_trans_L)'+sum(M_rot_L)'];
    FM_R = 1e-3*[sum(F_trans_R)'+sum(F_rot_R)'; sum(M_trans_R)'+sum(M_rot_R)'];
    
    FM_b = nan(6,1);
    
    FM_b(1:3)       = RL'*FM_L(1:3)+RR'*FM_R(1:3);
    FM_b(4:6)       = RL'*FM_L(4:6)+RR'*FM_R(4:6)+cross((Joint_left-cg_b),RL'*FM_L(1:3))+cross((Joint_right-cg_b),RR'*FM_R(1:3));
 
%     FM_strkpln(1:3) = R_strk*RL'*FM_L(1:3)+R_strk*RR'*FM_R(1:3);
%     FM_strkpln(4:6) = R_strk*RL'*FM_L(4:6)+R_strk*RR'*FM_R(4:6)+ ...
%                       R_strk*cross((Joint_left-cg_b),RL'*FM_L(1:3))+R_strk*cross((Joint_right-cg_b),RR'*FM_R(1:3));

    FM_strkpln(1:3) = R_strk*FM_b(1:3);
    FM_strkpln(4:6) = R_strk*FM_b(4:6);


end

