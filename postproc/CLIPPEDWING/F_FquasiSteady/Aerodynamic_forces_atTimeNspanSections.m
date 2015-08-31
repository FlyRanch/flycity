function [ FM_strkpln, FM_L, FM_R, U_left, U_right, Udot_left, Udot_right, alfa_L, alfa_R, alfa_dot_L, alfa_dot_R ] = Aerodynamic_forces_atTimeNspanSections( kine, body_model, wing_model, rot_lift_on)


    % Compute the aerodynamic forces on the wings by means of a
    % quasi-steady model:
    
    u_strk          = kine.u_strk;          % [ mm/s ]
    w_strk          = kine.w_strk;          % [ rad/s ]
    wL              = kine.wL;              % [ rad/s ]
    wR              = kine.wR;              % [ rad/s ]
    w_dot_L         = kine.w_dot_L;
    w_dot_R         = kine.w_dot_R;
    RL              = kine.RL;
    RR              = kine.RR;
        
    R_strk          = body_model.R_strk;
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
    
    N = size(u_strk,2);
    M = size(y_sect_L,2);
    
    nr_sect = length(y_sect_L);
    
    FM_strkpln  = zeros(6,N);
    FM_L        = zeros(6,N);
    FM_R        = zeros(6,N);
    
    FM_strkpln_trans  = zeros(6,N);
    FM_L_trans        = zeros(6,N);
    FM_R_trans        = zeros(6,N);
    
    FM_strkpln_rot    = zeros(6,N);
    FM_L_rot          = zeros(6,N);
    FM_R_rot          = zeros(6,N);
    
    UL          = zeros(3,nr_sect,N);                 % [ mm/s ]
    UR          = zeros(3,nr_sect,N);                 % [ mm/s ]
    U_dot_L     = zeros(3,nr_sect,N);                 % [ mm/s^2 ]
    U_dot_R     = zeros(3,nr_sect,N);                 % [ mm/s^2 ]
    alfa_L      = zeros(nr_sect,N);                   % [ rad ]
    alfa_R      = zeros(nr_sect,N);                   % [ rad ]
    alfa_dot_L  = zeros(nr_sect,N);                   % [ rad/s ]
    alfa_dot_R  = zeros(nr_sect,N);                   % [ rad/s ]
    L_trans_L   = zeros(nr_sect,N);  % [ N ]
    L_trans_R   = zeros(nr_sect,N);  % [ N ]
    D_trans_L   = zeros(nr_sect,N);  % [ N ]
    D_trans_R   = zeros(nr_sect,N);  % [ N ]
    F_trans_L   = zeros(nr_sect,3,N); % [ N ]
    F_trans_R   = zeros(nr_sect,3,N); % [ N ]
    M_trans_L   = zeros(nr_sect,3,N); % [ N*mm ]
    M_trans_R   = zeros(nr_sect,3,N); % [ N*mm ]
    F_rot_L     = zeros(nr_sect,3,N); % [ N ]
    F_rot_R     = zeros(nr_sect,3,N); % [ N ]
    M_rot_L     = zeros(nr_sect,3,N); % [ N*mm ]
    M_rot_R     = zeros(nr_sect,3,N); % [ N*mm ]
    
    vb = zeros(3,N);
    wb = zeros(3,N);


    for i = 1:N
        
        vb(:,i) = R_strk'*u_strk(:,i);
        wb(:,i) = R_strk'*w_strk(:,i);
        
        WL      = [ 0 -wL(3,i) wL(2,i); ...
                    wL(3,i) 0 -wL(1,i); ...
                   -wL(2,i) wL(1,i) 0];             % rad/s

        WR      = [ 0 -wR(3,i) wR(2,i); ...
                    wR(3,i) 0 -wR(1,i); ...
                   -wR(2,i) wR(1,i) 0];             % rad/s
    
        for j = 1:nr_sect
            
            % Compute the local left and right wing velocity at the 
            % current section:


            UL(:,j,i) = cross(wL(:,i),y_sect_L(:,j))+cross(RL(:,:,i)*wb(:,i),RL(:,:,i)*Joint_left+y_sect_L(:,j))+RL(:,:,i)*vb(:,i);   % [ mm/s ]
            UR(:,j,i) = cross(wR(:,i),y_sect_R(:,j))+cross(RR(:,:,i)*wb(:,i),RR(:,:,i)*Joint_right+y_sect_R(:,j))+RR(:,:,i)*vb(:,i);  % [ mm/s ]

            U_dot_L(:,j,i) = cross(w_dot_L(:,i),y_sect_L(:,j))+cross(WL*RL(:,:,i)*wb(:,i),RL(:,:,i)*Joint_left+y_sect_L(:,j))+...
                           cross(RL(:,:,i)*wb(:,i),WL*RL(:,:,i)*Joint_left+y_sect_L(:,j))+WL*RL(:,:,i)*vb(:,i);

            U_dot_R(:,j,i) = cross(w_dot_R(:,i),y_sect_R(:,j))+cross(WR*RR(:,:,i)*wb(:,i),RR(:,:,i)*Joint_right+y_sect_R(:,j))+...
                           cross(RR(:,:,i)*wb(:,i),WR*RR(:,:,i)*Joint_right+y_sect_R(:,j))+WR*RR(:,:,i)*vb(:,i);

            alfa_L(j,i) = real(atan2(UL(3,j,i),UL(1,j,i)));     % [ rad ]
            alfa_R(j,i) = real(atan2(UR(3,j,i),UR(1,j,i)));     % [ rad ]
%             alfa_dot_L(j,i) = -(UL(3,j,i)/(UL(1,j,i)^2+UL(3,j,i)^2))*U_dot_L(1)+(UL(1,j,i)/(UL(1,j,i)^2+UL(3,j,i)^2))*U_dot_L(3); % [ rad/s ]
%             alfa_dot_R(j,i) = -(UR(3,j,i)/(UR(1,j,i)^2+UR(3,j,i)^2))*U_dot_R(1)+(UR(1,j,i)/(UR(1,j,i)^2+UR(3,j,i)^2))*U_dot_R(3); % [ rad/s ]

            alfa_dot_L(j,i) = -wL(2,i);
            alfa_dot_R(j,i) = -wR(2,i);

            % Compute location center of pressure:
            
            x_cp_L = x_LE_L(j)-chords_L(j)*(0.82*abs(alfa_L(j,i))/pi+0.05);
            x_cp_R = x_LE_R(j)-chords_R(j)*(0.82*abs(alfa_R(j,i))/pi+0.05);

            if alfa_L(j,i) >= 0
                
                alfa_left        = abs(radtodeg(alfa_L(j,i)));                                                       % [ deg ]
                Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
                Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
                L_trans_L(j,i)   = 1e-3*0.5*rho*chords_L(j)*(UL(1,j,i)^2+UL(3,j,i)^2).*Cl_left*delta_R;         % [ N ]
                D_trans_L(j,i)   = 1e-3*0.5*rho*chords_L(j)*(UL(1,j,i)^2+UL(3,j,i)^2).*Cd_left*delta_R;         % [ N ]
                F_trans_L(j,:,i) = [(L_trans_L(j,i)*sind(alfa_left)-D_trans_L(j,i)*cosd(alfa_left)) ...
                                    0 ...
                                    -(L_trans_L(j,i)*cosd(alfa_left)+D_trans_L(j,i)*sind(alfa_left))];        % [ N ]
                M_trans_L(j,:,i) = (cross([x_cp_L; y_sect_L(2,j); y_sect_L(3,j)],F_trans_L(j,:,i)'))';                                   % [ N*mm ]
                
            else
                
                alfa_left        = abs(radtodeg(alfa_L(j,i)));                                                       % [ deg ]
                Cl_left          = 0.225+1.58*sind(2.13*alfa_left-7.2);
                Cd_left          = 1.92-1.55*cosd(2.04*alfa_left-9.82);
                L_trans_L(j,i)   = 1e-3*0.5*rho*chords_L(j)*(UL(1,j,i)^2+UL(3,j,i)^2).*Cl_left*delta_R;         % [ N ]
                D_trans_L(j,i)   = 1e-3*0.5*rho*chords_L(j)*(UL(1,j,i)^2+UL(3,j,i)^2).*Cd_left*delta_R;         % [ N ]
                F_trans_L(j,:,i) = [(L_trans_L(j,i).*sind(alfa_left)-D_trans_L(j,i).*cosd(alfa_left)) ...
                                    0 ...
                                    (L_trans_L(j,i).*cosd(alfa_left)+D_trans_L(j,i).*sind(alfa_left))];        % [ N ]
                M_trans_L(j,:,i) = (cross([x_cp_L; y_sect_L(2,j); y_sect_L(3,j)],F_trans_L(j,:,i)'))';                                   % [ N*mm ]
                
            end
            
            
            if alfa_R(j,i) >= 0
                
                alfa_right       = abs(radtodeg(alfa_R(j,i)));                                                       % [ deg ]
                Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
                Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
                L_trans_R(j,i)   = 1e-3*0.5*rho*chords_R(j)*(UR(1,j,i)^2+UR(3,j,i)^2).*Cl_right*delta_R;        % [ N ]
                D_trans_R(j,i)   = 1e-3*0.5*rho*chords_R(j)*(UR(1,j,i)^2+UR(3,j,i)^2).*Cd_right*delta_R;        % [ N ]
                F_trans_R(j,:,i) = [(L_trans_R(j,i)*sind(alfa_right)-D_trans_R(j,i)*cosd(alfa_right)) ...
                                    0 ...
                                    -(L_trans_R(j,i)*cosd(alfa_right)+D_trans_R(j,i)*sind(alfa_right))];      % [ N ]
                M_trans_R(j,:,i) = (cross([x_cp_R; y_sect_R(2,j); y_sect_R(3,j)],F_trans_R(j,:,i)'))';                                   % [ N*mm ]
                
            else
                
                alfa_right       = abs(radtodeg(alfa_R(j,i)));                                                       % [ deg ]
                Cl_right         = 0.225+1.58*sind(2.13*alfa_right-7.2);
                Cd_right         = 1.92-1.55*cosd(2.04*alfa_right-9.82);
                L_trans_R(j,i)   = 1e-3*0.5*rho*chords_R(j)*(UR(1,j,i)^2+UR(3,j,i)^2).*Cl_right*delta_R;        % [ N ]
                D_trans_R(j,i)   = 1e-3*0.5*rho*chords_R(j)*(UR(1,j,i)^2+UR(3,j,i)^2).*Cd_right*delta_R;        % [ N ]
                F_trans_R(j,:,i) = [(L_trans_R(j,i)*sind(alfa_right)-D_trans_R(j,i)*cosd(alfa_right)) ...
                                    0 ...
                                    (L_trans_R(j,i)*cosd(alfa_right)+D_trans_R(j,i)*sind(alfa_right))];      % [ N ]
                M_trans_R(j,:,i) = (cross([x_cp_R; y_sect_R(2,j); y_sect_R(3,j)],F_trans_R(j,:,i)'))';                                   % [ N*mm ]
                
            end
            
            if rot_lift_on == 1
            
                F_rot_L(j,:,i) = 1e-3*[ 0 0 C_rot*rho*alfa_dot_L(j,i).*chords_L(j)^2*sqrt(UL(1,j,i)^2+UL(3,j,i)^2)*delta_R];    % [ N ] 
                F_rot_R(j,:,i) = 1e-3*[ 0 0 C_rot*rho*alfa_dot_R(j,i).*chords_R(j)^2*sqrt(UR(1,j,i)^2+UR(3,j,i)^2)*delta_R];    % [ N ]

                M_rot_L(j,:,i) = (cross([x_cp_L; y_sect_L(2,j); y_sect_L(3,j)],F_rot_L(j,:,i)'))';                         % [ N*mm ]
                M_rot_R(j,:,i) = (cross([x_cp_R; y_sect_R(2,j); y_sect_R(3,j)],F_rot_R(j,:,i)'))';                         % [ N*mm ]    
            
            end
            
        end
    
    end
    

  %% sum values from span sections & Rstrk rotation
%      for i = 1:N
%             
%             FM_L(:,i) = [sum(F_trans_L(:,:,i))'+sum(F_rot_L(:,:,i))'; sum(M_trans_L(:,:,i))'+sum(M_rot_L(:,:,i))'];
%             FM_R(:,i) = [sum(F_trans_R(:,:,i))'+sum(F_rot_R(:,:,i))'; sum(M_trans_R(:,:,i))'+sum(M_rot_R(:,:,i))'];
%             
%             FM_L_trans(:,i)        = [sum(F_trans_L(:,:,i))'; sum(M_trans_L(:,:,i))'];
%             FM_R_trans(:,i)        = [sum(F_trans_R(:,:,i))'; sum(M_trans_R(:,:,i))'];
%             
%             FM_L_rot(:,i)          = [sum(F_rot_L(:,:,i))'; sum(M_rot_L(:,:,i))'];
%             FM_R_rot(:,i)          = [sum(F_rot_R(:,:,i))'; sum(M_rot_R(:,:,i))'];
%             
%             FM_strkpln(1:3,i) = R_strk*RL(:,:,i)'*FM_L(1:3,i)+R_strk*RR(:,:,i)'*FM_R(1:3,i);
%             FM_strkpln(4:6,i) = R_strk*RL(:,:,i)'*FM_L(4:6,i)+R_strk*RR(:,:,i)'*FM_R(4:6,i)+ ...
%                                 R_strk*cross((Joint_left-cg_b),(RL(:,:,i)'*FM_L(1:3,i)))+R_strk*cross((Joint_right-cg_b),(RR(:,:,i)'*FM_R(1:3,i)));
%             
%             FM_strkpln_trans(1:3,i) = R_strk*RL(:,:,i)'*FM_L_trans(1:3,i)+R_strk*RR(:,:,i)'*FM_R_trans(1:3,i);
%             FM_strkpln_trans(4:6,i) = R_strk*RL(:,:,i)'*FM_L_trans(4:6,i)+R_strk*RR(:,:,i)'*FM_R_trans(4:6,i)+ ...
%                                 R_strk*cross((Joint_left-cg_b),(RL(:,:,i)'*FM_L_trans(1:3,i)))+R_strk*cross((Joint_right-cg_b),(RR(:,:,i)'*FM_R_trans(1:3,i)));
%        
%             FM_strkpln_rot(1:3,i) = R_strk*RL(:,:,i)'*FM_L_rot(1:3,i)+R_strk*RR(:,:,i)'*FM_R_rot(1:3,i);
%             FM_strkpln_rot(4:6,i) = R_strk*RL(:,:,i)'*FM_L_rot(4:6,i)+R_strk*RR(:,:,i)'*FM_R_rot(4:6,i)+ ...
%                                 R_strk*cross((Joint_left-cg_b),(RL(:,:,i)'*FM_L_rot(1:3,i)))+R_strk*cross((Joint_right-cg_b),(RR(:,:,i)'*FM_R_rot(1:3,i)));
% 
% %             FM_strkpln(1:3,i) = R_strk*RL(:,:,i)'*FM_L(1:3,i)+R_strk*RR(:,:,i)'*FM_R(1:3,i);
% %             FM_strkpln(4:6,i) = R_strk*RL(:,:,i)'*FM_L(4:6,i)+R_strk*RR(:,:,i)'*FM_R(4:6,i)+ ...
% %                                 R_strk*cross([0; Joint_left(2); 0],(RL(:,:,i)'*FM_L(1:3,i)))+R_strk*cross([0; Joint_right(2); 0],(RR(:,:,i)'*FM_R(1:3,i)));
% %             
% %             FM_strkpln_trans(1:3,i) = R_strk*RL(:,:,i)'*FM_L_trans(1:3,i)+R_strk*RR(:,:,i)'*FM_R_trans(1:3,i);
% %             FM_strkpln_trans(4:6,i) = R_strk*RL(:,:,i)'*FM_L_trans(4:6,i)+R_strk*RR(:,:,i)'*FM_R_trans(4:6,i)+ ...
% %                                 R_strk*cross([0; Joint_left(2); 0],(RL(:,:,i)'*FM_L_trans(1:3,i)))+R_strk*cross([0; Joint_right(2); 0],(RR(:,:,i)'*FM_R_trans(1:3,i)));
% %        
% %             FM_strkpln_rot(1:3,i) = R_strk*RL(:,:,i)'*FM_L_rot(1:3,i)+R_strk*RR(:,:,i)'*FM_R_rot(1:3,i);
% %             FM_strkpln_rot(4:6,i) = R_strk*RL(:,:,i)'*FM_L_rot(4:6,i)+R_strk*RR(:,:,i)'*FM_R_rot(4:6,i)+ ...
% %                                 R_strk*cross([0; Joint_left(2); 0],(RL(:,:,i)'*FM_L_rot(1:3,i)))+R_strk*cross([0; Joint_right(2); 0],(RR(:,:,i)'*FM_R_rot(1:3,i)));
% % 
% %             FM_strkpln(1:3,i) = R_strk*RL(:,:,i)'*FM_L(1:3,i)+R_strk*RR(:,:,i)'*FM_R(1:3,i);
% %             FM_strkpln(4:6,i) = R_strk*RL(:,:,i)'*FM_L(4:6,i)+R_strk*RR(:,:,i)'*FM_R(4:6,i)+ ...
% %                                 R_strk*cross(Joint_left,(RL(:,:,i)'*FM_L(1:3,i)))+R_strk*cross(Joint_right,(RR(:,:,i)'*FM_R(1:3,i)));
% %             
% %             FM_strkpln_trans(1:3,i) = R_strk*RL(:,:,i)'*FM_L_trans(1:3,i)+R_strk*RR(:,:,i)'*FM_R_trans(1:3,i);
% %             FM_strkpln_trans(4:6,i) = R_strk*RL(:,:,i)'*FM_L_trans(4:6,i)+R_strk*RR(:,:,i)'*FM_R_trans(4:6,i)+ ...
% %                                 R_strk*cross(Joint_left,(RL(:,:,i)'*FM_L_trans(1:3,i)))+R_strk*cross(Joint_right,(RR(:,:,i)'*FM_R_trans(1:3,i)));
% %        
% %             FM_strkpln_rot(1:3,i) = R_strk*RL(:,:,i)'*FM_L_rot(1:3,i)+R_strk*RR(:,:,i)'*FM_R_rot(1:3,i);
% %             FM_strkpln_rot(4:6,i) = R_strk*RL(:,:,i)'*FM_L_rot(4:6,i)+R_strk*RR(:,:,i)'*FM_R_rot(4:6,i)+ ...
% %                                 R_strk*cross(Joint_left,(RL(:,:,i)'*FM_L_rot(1:3,i)))+R_strk*cross(Joint_right,(RR(:,:,i)'*FM_R_rot(1:3,i)));
% 
%      end
     
     
       %% force & torque rotation for all span sections (from wing coord system, to body coord system, to strokeplane coord system)
     for i = 1:N        % time loop
         for j = 1:M    % span section loop
             
            F_strkpln_trans_L(j,i,:) = R_strk*RL(:,:,i)'*F_trans_L(j,:,i)';
            M_strkpln_trans_L(j,i,:) = R_strk*RL(:,:,i)'*M_trans_L(j,:,i)';
            
            F_strkpln_rot_L(j,i,:) = R_strk*RL(:,:,i)'*F_rot_L(j,:,i)';
            M_strkpln_rot_L(j,i,:) = R_strk*RL(:,:,i)'*M_rot_L(j,:,i)';
            
            F_strkpln_trans_R(j,i,:) = R_strk*RR(:,:,i)'*F_trans_R(j,:,i)';
            M_strkpln_trans_R(j,i,:) = R_strk*RR(:,:,i)'*M_trans_R(j,:,i)';
            
            F_strkpln_rot_R(j,i,:) = R_strk*RR(:,:,i)'*F_rot_R(j,:,i)';
            M_strkpln_rot_R(j,i,:) = R_strk*RR(:,:,i)'*M_rot_R(j,:,i)';
         end
     end
     
     
     
    %% switch y,z position in velocity array
     U_left = zeros(3,N,nr_sect);
     U_right = zeros(3,N,nr_sect);
     
     Udot_left = zeros(3,N,nr_sect);
     Udot_right = zeros(3,N,nr_sect);
     
     for i = 1:N
         
         for j = 1:nr_sect
         
            U_left(:,i,j) = UL(:,j,i);
            U_right(:,i,j) = UR(:,j,i);
         
            Udot_left(:,i,j) = U_dot_L(:,j,i);
            Udot_right(:,i,j) = U_dot_R(:,j,i);
         
         end
         
     end    

     
     %% store data
     
     %% strokeplane coord system
     % translational forces
     FM_L.Fx_strkpln_trans_L(:,:) = F_strkpln_trans_L(:,:,1);
     FM_L.Fy_strkpln_trans_L(:,:) = F_strkpln_trans_L(:,:,2);
     FM_L.Fz_strkpln_trans_L(:,:) = F_strkpln_trans_L(:,:,3);
     
     FM_L.Mx_strkpln_trans_L(:,:) = M_strkpln_trans_L(:,:,1);
     FM_L.My_strkpln_trans_L(:,:) = M_strkpln_trans_L(:,:,2);
     FM_L.Mz_strkpln_trans_L(:,:) = M_strkpln_trans_L(:,:,3);
     
     FM_R.Fx_strkpln_trans_R(:,:) = F_strkpln_trans_R(:,:,1);
     FM_R.Fy_strkpln_trans_R(:,:) = F_strkpln_trans_R(:,:,2);
     FM_R.Fz_strkpln_trans_R(:,:) = F_strkpln_trans_R(:,:,3);
     
     FM_R.Mx_strkpln_trans_R(:,:) = M_strkpln_trans_R(:,:,1);
     FM_R.My_strkpln_trans_R(:,:) = M_strkpln_trans_R(:,:,2);
     FM_R.Mz_strkpln_trans_R(:,:) = M_strkpln_trans_R(:,:,3);
     
     % rotational forces
     FM_L.Fx_strkpln_rot_L(:,:) = F_strkpln_rot_L(:,:,1);
     FM_L.Fy_strkpln_rot_L(:,:) = F_strkpln_rot_L(:,:,2);
     FM_L.Fz_strkpln_rot_L(:,:) = F_strkpln_rot_L(:,:,3);
     
     FM_L.Mx_strkpln_rot_L(:,:) = M_strkpln_rot_L(:,:,1);
     FM_L.My_strkpln_rot_L(:,:) = M_strkpln_rot_L(:,:,2);
     FM_L.Mz_strkpln_rot_L(:,:) = M_strkpln_rot_L(:,:,3);
     
     FM_R.Fx_strkpln_rot_R(:,:) = F_strkpln_rot_R(:,:,1);
     FM_R.Fy_strkpln_rot_R(:,:) = F_strkpln_rot_R(:,:,2);
     FM_R.Fz_strkpln_rot_R(:,:) = F_strkpln_rot_R(:,:,3);
     
     FM_R.Mx_strkpln_rot_R(:,:) = M_strkpln_rot_R(:,:,1);
     FM_R.My_strkpln_rot_R(:,:) = M_strkpln_rot_R(:,:,2);
     FM_R.Mz_strkpln_rot_R(:,:) = M_strkpln_rot_R(:,:,3);

     % translational+rotational forces
     FM_L.Fx_strkpln_L(:,:) = F_strkpln_trans_L(:,:,1) +  F_strkpln_rot_L(:,:,1);
     FM_L.Fy_strkpln_L(:,:) = F_strkpln_trans_L(:,:,2) +  F_strkpln_rot_L(:,:,2);
     FM_L.Fz_strkpln_L(:,:) = F_strkpln_trans_L(:,:,3) +  F_strkpln_rot_L(:,:,3);
     
     FM_R.Fx_strkpln_R(:,:) = F_strkpln_trans_R(:,:,1) +  F_strkpln_rot_R(:,:,1);
     FM_R.Fy_strkpln_R(:,:) = F_strkpln_trans_R(:,:,2) +  F_strkpln_rot_R(:,:,2);
     FM_R.Fz_strkpln_R(:,:) = F_strkpln_trans_R(:,:,3) +  F_strkpln_rot_R(:,:,3);
     
     FM_L.Mx_strkpln_L(:,:) = M_strkpln_trans_L(:,:,1) +  M_strkpln_rot_L(:,:,1);
     FM_L.My_strkpln_L(:,:) = M_strkpln_trans_L(:,:,2) +  M_strkpln_rot_L(:,:,2);
     FM_L.Mz_strkpln_L(:,:) = M_strkpln_trans_L(:,:,3) +  M_strkpln_rot_L(:,:,3);
     
     FM_R.Mx_strkpln_R(:,:) = M_strkpln_trans_R(:,:,1) +  M_strkpln_rot_R(:,:,1);
     FM_R.My_strkpln_R(:,:) = M_strkpln_trans_R(:,:,2) +  M_strkpln_rot_R(:,:,2);
     FM_R.Mz_strkpln_R(:,:) = M_strkpln_trans_R(:,:,3) +  M_strkpln_rot_R(:,:,3);
     
     % mean forces & torques
     F_mean(:,:) = [sum(F_strkpln_trans_L)+sum(F_strkpln_rot_L)+sum(F_strkpln_trans_R)+sum(F_strkpln_rot_R)];
     M_mean(:,:) = [sum(M_strkpln_trans_L)+sum(M_strkpln_rot_L)+sum(M_strkpln_trans_R)+sum(M_strkpln_rot_R)];
     FM_strkpln = [F_mean M_mean];

     %% wing coord system
     % translational forces
     FM_L.Fx_trans_L(:,:) = F_trans_L(:,1,:);
     FM_L.Fy_trans_L(:,:) = F_trans_L(:,2,:);
     FM_L.Fz_trans_L(:,:) = F_trans_L(:,3,:);
     
     FM_L.Mx_trans_L(:,:) = M_trans_L(:,1,:);
     FM_L.My_trans_L(:,:) = M_trans_L(:,2,:);
     FM_L.Mz_trans_L(:,:) = M_trans_L(:,3,:);
     
     FM_R.Fx_trans_R(:,:) = F_trans_R(:,1,:);
     FM_R.Fy_trans_R(:,:) = F_trans_R(:,2,:);
     FM_R.Fz_trans_R(:,:) = F_trans_R(:,3,:);
     
     FM_R.Mx_trans_R(:,:) = M_trans_R(:,1,:);
     FM_R.My_trans_R(:,:) = M_trans_R(:,2,:);
     FM_R.Mz_trans_R(:,:) = M_trans_R(:,3,:);
     
     % rotational forces
     FM_L.Fx_rot_L(:,:) = F_rot_L(:,1,:);
     FM_L.Fy_rot_L(:,:) = F_rot_L(:,2,:);
     FM_L.Fz_rot_L(:,:) = F_rot_L(:,3,:);
     
     FM_L.Mx_rot_L(:,:) = M_rot_L(:,1,:);
     FM_L.My_rot_L(:,:) = M_rot_L(:,2,:);
     FM_L.Mz_rot_L(:,:) = M_rot_L(:,3,:);
     
     FM_R.Fx_rot_R(:,:) = F_rot_R(:,1,:);
     FM_R.Fy_rot_R(:,:) = F_rot_R(:,2,:);
     FM_R.Fz_rot_R(:,:) = F_rot_R(:,3,:);
     
     FM_R.Mx_rot_R(:,:) = M_rot_R(:,1,:);
     FM_R.My_rot_R(:,:) = M_rot_R(:,2,:);
     FM_R.Mz_rot_R(:,:) = M_rot_R(:,3,:);

     % translational+rotational forces
     FM_L.Fx_L(:,:) = F_trans_L(:,1,:) +  F_rot_L(:,1,:);
     FM_L.Fy_L(:,:) = F_trans_L(:,2,:) +  F_rot_L(:,2,:);
     FM_L.Fz_L(:,:) = F_trans_L(:,3,:) +  F_rot_L(:,3,:);
     
     FM_R.Fx_R(:,:) = F_trans_R(:,1,:) +  F_rot_R(:,1,:);
     FM_R.Fy_R(:,:) = F_trans_R(:,2,:) +  F_rot_R(:,2,:);
     FM_R.Fz_R(:,:) = F_trans_R(:,3,:) +  F_rot_R(:,3,:);
     
     FM_L.Mx_L(:,:) = M_trans_L(:,1,:) +  M_rot_L(:,1,:);
     FM_L.My_L(:,:) = M_trans_L(:,2,:) +  M_rot_L(:,2,:);
     FM_L.Mz_L(:,:) = M_trans_L(:,3,:) +  M_rot_L(:,3,:);
     
     FM_R.Mx_R(:,:) = M_trans_R(:,1,:) +  M_rot_R(:,1,:);
     FM_R.My_R(:,:) = M_trans_R(:,2,:) +  M_rot_R(:,2,:);
     FM_R.Mz_R(:,:) = M_trans_R(:,3,:) +  M_rot_R(:,3,:);
     
     % mean forces & torques
     F_L_mean(:,:) = sum(F_trans_L)+sum(F_rot_L);
     M_L_mean(:,:) = sum(M_trans_L)+sum(M_rot_L);
     FM_L.FM_mean = [F_L_mean M_L_mean];

     F_R_mean(:,:) = sum(F_trans_R)+sum(F_rot_R);
     M_R_mean(:,:) = sum(M_trans_R)+sum(M_rot_R);
     FM_R.FM_mean = [F_R_mean M_R_mean];
end

