function [ FM_strkpln, FM_L, FM_R ,U_left, U_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aerodynamic_forces( kine, body_model, wing_model, wb, rot_lift_on)


    % Compute the aerodynamic forces on the wings by means of a
    % quasi-steady model:
    
    u_strk          = kine.u_strk;          % [ mm/s ]
    w_strk          = kine.w_strk;          % [ rad/s ]
    wL              = kine.wL;              % [ rad/s ]
    wR              = kine.wR;              % [ rad/s ]
    RL              = kine.RL;
    RR              = kine.RR;
    R_strk          = kine.R_strk;
    
    Joint_left      = body_model.Joint_left;    % [ mm ]
    Joint_right     = body_model.Joint_right;   % [ mm ]
    cg_b            = body_model.cg_b;          % [ mm ]
    
    y_sect_L        = wing_model.y_sect_L;      % [ mm ]
    chords_L        = wing_model.chords_L;      % [ mm ]
    y_sect_R        = wing_model.y_sect_R;      % [ mm ]
    chords_R        = wing_model.chords_R;      % [ mm ]
    rho             = wing_model.rho;           % [ kg/mm^3 ]
    
    wb_loc          = wb.wb_loc;
    down_loc        = wb.down_loc;
    up_loc          = wb.up_loc;
    dt              = wb.dt;
    
    delta_R         = abs(y_sect_L(2,2)-y_sect_L(2,1));     % [ mm ]
    
    C_rot           = 1.55;
    
    N = size(u_strk,2);
    
    nr_sect = length(y_sect_L);
    
    FM_strkpln  = nan(6,N);
    FM_L        = nan(6,N);
    FM_R        = nan(6,N);
    
    UL          = nan(3,nr_sect,N);                 % [ mm/s ]
    UR          = nan(3,nr_sect,N);                 % [ mm/s ]
    alfa_L      = nan(nr_sect,N);                   % [ rad ]
    alfa_R      = nan(nr_sect,N);                   % [ rad ]
    alfa_dot_L  = nan(nr_sect,N);                   % [ rad/s ]
    alfa_dot_R  = nan(nr_sect,N);                   % [ rad/s ]
    
    for i = 1:N
            
            for j = 1:nr_sect
                
                % Compute the local left and right wing velocity at the 
                % current section:
                
                UL(:,j,i) = cross(wL(:,i),y_sect_L(:,j))+cross(RL(:,:,i)*R_strk'*w_strk(:,i),RL(:,:,i)*Joint_left+y_sect_L(:,j))+RL(:,:,i)*R_strk'*u_strk(:,i);   % [ mm/s ]
                UR(:,j,i) = cross(wR(:,i),y_sect_R(:,j))+cross(RR(:,:,i)*R_strk'*w_strk(:,i),RR(:,:,i)*Joint_right+y_sect_R(:,j))+RR(:,:,i)*R_strk'*u_strk(:,i);  % [ mm/s ]
                
                alfa_L(j,i) = real(atan2(UL(3,j,i),UL(1,j,i)));     % [ rad ]
                alfa_R(j,i) = real(atan2(UR(3,j,i),UR(1,j,i)));     % [ rad ]
                
            end
            
    end
    
    alfa_dot_L_t = nan(nr_sect,N-1);
    alfa_dot_R_t = nan(nr_sect,N-1);
    
    for k = 1:nr_sect
    
        alfa_dot_L_t(k,:) = diff(alfa_L(k,:));
        
        alfa_dot_R_t(k,:) = diff(alfa_R(k,:));
    
    end
    
    for i = 1:N
        
        for j = 1:nr_sect
            
            if i == 1
                
                if abs(alfa_dot_L_t(j,i)) >= (pi/2)
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i+1)/dt;
                    
                else
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i)/dt;
                    
                end
                
                if abs(alfa_dot_R_t(j,i)) >= (pi/2)
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i+1)/dt;
                    
                else
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i)/dt;
                    
                end
                
            elseif i == 2
                
                if abs(alfa_dot_L_t(j,i)) >= (pi/2)
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i+1)/dt;
                    
                else
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i)/dt;
                    
                end
                
                if abs(alfa_dot_R_t(j,i)) >= (pi/2)
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i+1)/dt;
                    
                else
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i)/dt;
                    
                end
                
            elseif i == (N-1)
                
                if abs(alfa_dot_L_t(j,i)) >= (pi/2)
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i-1)/dt;
                    
                else
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i)/dt;
                    
                end
                
                if abs(alfa_dot_R_t(j,i)) >= (pi/2)
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i-1)/dt;
                    
                else
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i)/dt;
                    
                end
                
            elseif i == N
                
                if abs(alfa_dot_L_t(j,i-1)) >= (pi/2)
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i-2)/dt;
                    
                else
                    
                    alfa_dot_L(j,i) = alfa_dot_L_t(j,i-1)/dt;
                    
                end
                
                if abs(alfa_dot_R_t(j,i-1)) >= (pi/2)
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i-2)/dt;
                    
                else
                    
                    alfa_dot_R(j,i) = alfa_dot_R_t(j,i-1)/dt;
                    
                end
               
            elseif i > 2 && i < (N-1)
            
                if abs(alfa_dot_L_t(j,i)) >= (pi/2)

                    alfa_dot_L(j,i) = (2/3)*alfa_dot_L_t(j,i-1)/dt+(1/3)*alfa_dot_L_t(j,i+1)/dt;

                elseif abs(alfa_dot_L_t(j,i-1)) >= (pi/2)

                    alfa_dot_L(j,i) = (1/3)*alfa_dot_L_t(j,i-2)/dt+(2/3)*alfa_dot_L_t(j,i)/dt;

                else
                    
                    alfa_dot_L(j,i) = 0.5*alfa_dot_L_t(j,i-1)/dt+0.5*alfa_dot_L_t(j,i)/dt;
                    
                end


                if abs(alfa_dot_R_t(j,i)) >= (pi/2)

                    alfa_dot_R(j,i) = (2/3)*alfa_dot_R_t(j,i-1)/dt+(1/3)*alfa_dot_R_t(j,i+1)/dt;

                elseif abs(alfa_dot_R_t(j,i-1)) >= (pi/2)

                    alfa_dot_R(j,i) = (1/3)*alfa_dot_R_t(j,i-2)/dt+(2/3)*alfa_dot_R_t(j,i)/dt;

                else
                    
                    alfa_dot_R(j,i) = 0.5*alfa_dot_R_t(j,i-1)/dt+0.5*alfa_dot_R_t(j,i)/dt;
                    
                end
             
            end
        
        end
    
    
    end
    
    L_trans_L  = zeros(nr_sect,N);  % [ N ]
    L_trans_R  = zeros(nr_sect,N);  % [ N ]
    D_trans_L  = zeros(nr_sect,N);  % [ N ]
    D_trans_R  = zeros(nr_sect,N);  % [ N ]
    F_trans_L = zeros(nr_sect,3,N); % [ N ]
    F_trans_R = zeros(nr_sect,3,N); % [ N ]
    M_trans_L = zeros(nr_sect,3,N); % [ N*mm ]
    M_trans_R = zeros(nr_sect,3,N); % [ N*mm ]
    F_rot_L   = zeros(nr_sect,3,N); % [ N ]
    F_rot_R   = zeros(nr_sect,3,N); % [ N ]
    M_rot_L   = zeros(nr_sect,3,N); % [ N*mm ]
    M_rot_R   = zeros(nr_sect,3,N); % [ N*mm ]
    
    for i = 1:N
            
            if i >= down_loc(1) && i <= down_loc(2)
                
                alfa_left        = radtodeg(alfa_L(:,i));                                                       % [ deg ]
                alfa_right       = radtodeg(alfa_R(:,i));                                                       % [ deg ]
                Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
                Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
                Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
                Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
                L_trans_L(:,i)   = 1e-3*0.5*rho*chords_L.*((UL(1,:,i)').^2+(UL(3,:,i)').^2).*Cl_left*delta_R;        % [ N ]
                L_trans_R(:,i)   = 1e-3*0.5*rho*chords_R.*((UR(1,:,i)').^2+(UR(3,:,i)').^2).*Cl_right*delta_R;       % [ N ]
                D_trans_L(:,i)   = 1e-3*0.5*rho*chords_L.*((UL(1,:,i)').^2+(UL(3,:,i)').^2).*Cd_left*delta_R;        % [ N ]
                D_trans_R(:,i)   = 1e-3*0.5*rho*chords_R.*((UR(1,:,i)').^2+(UR(3,:,i)').^2).*Cd_right*delta_R;       % [ N ]
                F_trans_L(:,:,i) = [(L_trans_L(:,i).*sind(alfa_left)-D_trans_L(:,i).*cosd(alfa_left)) ...
                                    zeros(nr_sect,1) ...
                                    -(L_trans_L(:,i).*cosd(alfa_left)+D_trans_L(:,i).*sind(alfa_left))];        % [ N ]
                F_trans_R(:,:,i) = [(L_trans_R(:,i).*sind(alfa_right)-D_trans_R(:,i).*cosd(alfa_right)) ...
                                    zeros(nr_sect,1) ...
                                    -(L_trans_R(:,i).*cosd(alfa_right)+D_trans_R(:,i).*sind(alfa_right))];      % [ N ]
                
                    for j = 1:nr_sect

                        M_trans_L(j,:,i) = (cross(y_sect_L(:,j),F_trans_L(j,:,i)'))';                     % [ N*mm ]
                        M_trans_R(j,:,i) = (cross(y_sect_R(:,j),F_trans_R(j,:,i)'))';                     % [ N*mm ]

                    end
                
                if rot_lift_on == 1
                    
                    F_rot_L(:,:,i) = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_L(:,i).*chords_L.^2.*sqrt((UL(1,:,i)').^2+(UL(3,:,i)').^2)*delta_R];     % [ N ]
                    F_rot_R(:,:,i) = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_R(:,i).*chords_R.^2.*sqrt((UR(1,:,i)').^2+(UR(3,:,i)').^2)*delta_R];     % [ N ]
                    
                    for j = 1:nr_sect

                        M_rot_L(j,:,i) = (cross(y_sect_L(:,j),F_rot_L(j,:,i)'))';                         % [ N*mm ]
                        M_rot_R(j,:,i) = (cross(y_sect_R(:,j),F_rot_R(j,:,i)'))';                         % [ N*mm ]

                    end
                    
                end
                
            elseif i >= up_loc(1) && i <= up_loc(2)
                
                alfa_left        = -radtodeg(alfa_L(:,i));                                                      % [ N ] 
                alfa_right       = -radtodeg(alfa_R(:,i));                                                      % [ N ] 
                Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
                Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
                Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
                Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
                L_trans_L(:,i)   = 1e-3*0.5*rho*chords_L.*((UL(1,:,i)').^2+(UL(3,:,i)').^2).*Cl_left*delta_R;        % [ N ] 
                L_trans_R(:,i)   = 1e-3*0.5*rho*chords_R.*((UR(1,:,i)').^2+(UR(3,:,i)').^2).*Cl_right*delta_R;       % [ N ] 
                D_trans_L(:,i)   = 1e-3*0.5*rho*chords_L.*((UL(1,:,i)').^2+(UL(3,:,i)').^2).*Cd_left*delta_R;        % [ N ] 
                D_trans_R(:,i)   = 1e-3*0.5*rho*chords_R.*((UR(1,:,i)').^2+(UR(3,:,i)').^2).*Cd_right*delta_R;       % [ N ] 
                F_trans_L(:,:,i)  = [L_trans_L(:,i).*sind(alfa_left)-D_trans_L(:,i).*cosd(alfa_left) ...
                                     zeros(nr_sect,1) ...
                                     L_trans_L(:,i).*cosd(alfa_left)+D_trans_L(:,i).*sind(alfa_left)];          % [ N ] 
                F_trans_R(:,:,i)  = [L_trans_R(:,i).*sind(alfa_right)-D_trans_R(:,i).*cosd(alfa_right) ...
                                     zeros(nr_sect,1) ...
                                     L_trans_R(:,i).*cosd(alfa_right)+D_trans_R(:,i).*sind(alfa_right)];        % [ N ] 
                                 
                for j = 1:nr_sect

                        M_trans_L(j,:,i) = (cross(y_sect_L(:,j),F_trans_L(j,:,i)'))';                     % [ N*mm ]
                        M_trans_R(j,:,i) = (cross(y_sect_R(:,j),F_trans_R(j,:,i)'))';                     % [ N*mm ]

                end
                
                if rot_lift_on == 1
                    
                    F_rot_L(:,:,i) = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_L(:,i).*chords_L.^2.*sqrt((UL(1,:,i)').^2+(UL(3,:,i)').^2)*delta_R];    % [ N ] 
                    F_rot_R(:,:,i) = 1e-3*[ zeros(nr_sect,1) zeros(nr_sect,1) C_rot*rho*alfa_dot_R(:,i).*chords_R.^2.*sqrt((UR(1,:,i)').^2+(UR(3,:,i)').^2)*delta_R];    % [ N ]
                    
                    for j = 1:nr_sect

                        M_rot_L(j,:,i) = (cross(y_sect_L(:,j),F_rot_L(j,:,i)'))';                         % [ N*mm ]
                        M_rot_R(j,:,i) = (cross(y_sect_R(:,j),F_rot_R(j,:,i)'))';                         % [ N*mm ]

                    end
                    
                end
                
                
                
            end
            
    end
    
%     UL_t = nan(nr_sect,N);
%     UR_t = nan(nr_sect,N);
%     
%     UL_t(:,1:N) = UL(3,:,1:N);
%     UR_t(:,1:N) = UR(3,:,1:N);
%     
%     figure()
%     surf(UL_t)
%     
%     figure()
%     surf(UR_t)
%     
%     figure()
%     surf(alfa_L)
%     
%     figure()
%     surf(alfa_R)
%     
%     figure()
%     surf(L_trans_L)
%     
%     figure()
%     surf(L_trans_R)
%     
%     figure()
%     surf(D_trans_L)
%     
%     figure()
%     surf(D_trans_R)
%     
%     F_t_L = nan(nr_sect,N);
%     F_t_R = nan(nr_sect,N);
%     F_r_L = nan(nr_sect,N);
%     F_r_R = nan(nr_sect,N);
%     
%     F_t_L(:,1:N) = F_trans_L(:,3,1:N);
%     F_t_R(:,1:N) = F_trans_R(:,3,1:N);
%     F_r_L(:,1:N) = F_rot_L(:,3,1:N);
%     F_r_R(:,1:N) = F_rot_R(:,3,1:N);
%     
%     figure()
%     surf(F_t_L)
%     
%     figure()
%     surf(F_t_R)
%     
%     figure()
%     surf(F_r_L)
%     
%     figure()
%     surf(F_r_R)
    
     for i = 1:N
            
            FM_L(:,i) = [sum(F_trans_L(:,:,i))'+sum(F_rot_L(:,:,i))'; sum(M_trans_L(:,:,i))'+sum(M_rot_L(:,:,i))'];
            FM_R(:,i) = [sum(F_trans_R(:,:,i))'+sum(F_rot_R(:,:,i))'; sum(M_trans_R(:,:,i))'+sum(M_rot_R(:,:,i))'];
            
            FM_strkpln(1:3,i) = R_strk*RL(:,:,i)'*FM_L(1:3,i)+R_strk*RR(:,:,i)'*FM_R(1:3,i);
            FM_strkpln(4:6,i) = R_strk*RL(:,:,i)'*FM_L(4:6,i)+R_strk*RR(:,:,i)'*FM_R(4:6,i); %+ ...
                                %R_strk*cross((Joint_left-cg_b),RL(:,:,i)'*FM_L(1:3,i))+R_strk*cross((Joint_right-cg_b),RR(:,:,i)'*FM_R(1:3,i));
            
       
     end
    
     U_left = nan(3,N,nr_sect);
     U_right = nan(3,N,nr_sect);
     
     for i = 1:N
         
         for j = 1:nr_sect
         
            U_left(:,i,j) = UL(:,j,i);
            U_right(:,i,j) = UR(:,j,i);
         
         end
         
     end
    

end

