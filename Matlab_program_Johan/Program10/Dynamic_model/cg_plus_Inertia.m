function [cg_body, cg_L, cg_R, I_body, I_wing, I_v_wing, m_b, m_w, m_v_w ,y_sect, chords] = cg_plus_Inertia(settings, pathDB, m_wl, w_l ,seq_nr, nr_sect )


    % Load body model:

    [x, y, z] = load_body_model(settings, seq_nr, [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1] );
    
    
    
    % Set number of elements, position of the elements and in x and y direction for body and wing
    
    sect_nr_b_x = 100;
    
    sect_nr_b_y = 50;
    
    sect_nr_w_x = 100;
    
    sect_nr_w_y = 100;
    
    
    temp_x_cont_b = x{1};
    
    temp_y_cont_b = y{1};
    
    temp_x_cont_w = x{3};
    
    temp_y_cont_w = y{3};
    
    
    J_pos_L = pathDB.joint_pos_L(:,seq_nr);
    
    J_pos_R = pathDB.joint_pos_R(:,seq_nr);
    
    
    x_cont_b = [temp_x_cont_b(1:15,7); temp_x_cont_b(17:21,7)];
    
    y_cont_b = [temp_y_cont_b(1:15,7); temp_y_cont_b(17:21,7)];
    
    x_cont_w = -J_pos_R(1) + temp_x_cont_w(:,1);
    
    y_cont_w = -J_pos_R(2) + temp_y_cont_w(:,1);
    
    clear temp_x_cont_b temp_y_cont_b temp_x_cont_w temp_y_cont_w
    
    
    
    b_length = pathDB.body_l(seq_nr);
    
    b_width = 2*max(y_cont_b);
    
    
    w_length = pathDB.wing_l(seq_nr);
    
    w_width = max(x_cont_w)-min(x_cont_w);
    
    
    delta_x_b = b_length/sect_nr_b_x;
    
    delta_y_b = b_width/sect_nr_b_y;
    
    delta_x_w = w_width/sect_nr_w_x;
    
    delta_y_w = w_length/sect_nr_w_y;
    
    
    pos_b_x = min(x_cont_b)+((0.5/(sect_nr_b_x)):(1/(sect_nr_b_x)):(1-(0.5/(sect_nr_b_x)))).*b_length;
    
    pos_b_y = ((0.5/(sect_nr_b_y)):(1/(sect_nr_b_y)):(1-(0.5/(sect_nr_b_y)))).*b_width;
          
    pos_w_x = min(x_cont_w)+((0.5/(sect_nr_w_x)):(1/(sect_nr_w_x)):(1-(0.5/(sect_nr_w_x)))).*w_width;
    
    pos_w_y = ((0.5/(sect_nr_w_y)):(1/(sect_nr_w_y)):(1-(0.5/(sect_nr_w_y)))).*w_length;
    

    % Interpolate coordinates on the contours of the wing and the body:
    
    y_sect_b = interp1(x_cont_b, y_cont_b, pos_b_x);
    
    x_sect_1_w = interp1(y_cont_w(1:((end+1)/2)),x_cont_w(1:((end+1)/2)),pos_w_y);
    
    x_sect_2_w = interp1(y_cont_w(((end+1)/2):end),x_cont_w(((end+1)/2):end),pos_w_y);
    
    [x_min, N_x_min] = min(x_cont_w);
    
    [x_max, N_x_max] = max(x_cont_w);
    
    N_x = length(x_cont_w);
    
    y_sect_1_w = interp1(x_cont_w(N_x_max:N_x_min),y_cont_w(N_x_max:N_x_min),pos_w_x);
    
    y_sect_2_w = interp1(x_cont_w([N_x:-1:N_x_min 1:N_x_max]),y_cont_w([N_x:-1:N_x_min 1:N_x_max]),pos_w_x);
    
    clear x_min x_max N_x_Min N_x_max
    
    
    % Determine volume body:
    
    vol_b = sum(2*pi*0.5*y_sect_b.*y_sect_b*delta_x_b); % [mm^3]
    
    
    % Determine area wing:
    
    area_w = sum((x_sect_1_w-x_sect_2_w)*delta_y_w); % [mm^2]
    
    
    % Determine virtual air volume wing:
    
    vol_v_w = sum(pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w);
    
   
    % Determine air density:
    
    pres = 101325; % Pressure in Pascal
    
    Temp = 297.039; % Average temperature in Kelvin
    
    R = 287.058; % Gas constant for air
    
    rho_air = 1e-9*pres/(R*Temp); % Density in [kg/mm^3]
    
    
    % Calculate wing and body mass:
       
    m_w = w_l*area_w*1e-9*1200*5.4e-4;   % [kg]
    
    m_b = m_wl - 2*m_w;   % [kg]
    
    m_v_w = vol_v_w*rho_air;   % [kg]
    
   
    % Determine body density:
    
    rho_body = m_b/vol_b; % density in [kg/mm^3]
    

    % Determine wing density (2D):
    
    rho_w = m_w/area_w; % density (2D) in [kg/mm^2]
    
    
    % Determine the center of gravity of the body:
    
%     x_b_cg = sum((y_sect_b.*delta_x_b).*pos_b_x)/sum(y_sect_b.*delta_x_b);

    x_b_cg = 0; % Set to zero as the body model misses the abodomen often.
    
    % Determine the center of gravity of the wing:
    
    x_w_cg = sum(((y_sect_1_w-y_sect_2_w)*delta_x_w).*pos_w_x)/area_w;
    y_w_cg = sum(((x_sect_1_w-x_sect_2_w)*delta_y_w).*pos_w_y)/area_w;
    
    
    % Determine the center of gravity of the wing including virtual wing
    % mass:
    
    % Neglectable: effect of virtual mass is at most a 1.2% shift in the x
    % position and a .6 % change in y position
    
%     x_w_v_cg = sum((x_sect_1_w-0.5*(x_sect_1_w-x_sect_2_w)).*(rho_air*pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w)+(y_sect_1_w-y_sect_2_w)*delta_x_w*rho_w.*pos_w_x)/sum(rho_air*pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w+(x_sect_1_w-x_sect_2_w)*delta_y_w*rho_w);
%     y_w_v_cg = sum((rho_air*pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w+(x_sect_1_w-x_sect_2_w)*delta_y_w*rho_w).*pos_w_y)/sum(rho_air*pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w+(x_sect_1_w-x_sect_2_w)*delta_y_w*rho_w);

    
    % Determine the inertia tensor of the body:
    
    Ixx_b = sum(0.5*pi*rho_body*delta_x_b*y_sect_b.^4); % [kg * mm^2]
    
    Iyy_b = sum((1/12)*pi*rho_body*delta_x_b*y_sect_b.^2.*(3*y_sect_b.^2+delta_x_b^2)+pi*rho_body*delta_x_b*y_sect_b.^2.*abs(pos_b_x-x_b_cg).^2); % [kg * mm^2]
    
    Izz_b = Iyy_b; % [kg * mm^2]
    
    I_body = [ Ixx_b 0 0; ...
               0 Iyy_b 0; ...
               0 0 Izz_b];
           
    
    % Determine the inertia tensor of the wing:
    
    Ixx_w = sum(delta_y_w*rho_w*(x_sect_1_w-x_sect_2_w).*(pos_w_y - y_w_cg).^2); % [kg * mm^2]
    
    Iyy_w = sum(delta_x_w*rho_w*(y_sect_1_w-y_sect_2_w).*(pos_w_x - x_w_cg).^2); % [kg * mm^2]
    
    Izz_w = sum((1/12)*delta_y_w*rho_w*(x_sect_1_w-x_sect_2_w).*((x_sect_1_w-x_sect_2_w).^2+delta_y_w.^2)+delta_y_w*rho_w*(x_sect_1_w-x_sect_2_w).*((pos_w_y - y_w_cg).^2+(pos_w_x - x_w_cg).^2)); % [kg * mm^2]
    
    %Ixy_w = sum(-0.25*rho_w*((x_sect_1_w-x_w_cg).^2-(x_sect_2_w-x_w_cg).^2).*((pos_w_y+0.5*delta_y_w-y_w_cg).^2-(pos_w_y-0.5*delta_y_w-y_w_cg).^2)); % [kg * mm^2]
    
    %Ixy_w = -sum((rho_w*delta_y_w*(x_sect_1_w-x_sect_2_w).*(pos_w_y-y_w_cg)).*(rho_w*delta_x_w*(y_sect_1_w-y_sect_2_w).*(pos_w_x-x_w_cg)))
    
    I_wing = [ Ixx_w 0 0; ...
               0 Iyy_w 0; ...
               0 0 Izz_w ]; % [ kg * mm^2 ]
           
           
    % Determine the inertia tensor of the wing including virtual wing mass:
    
    Ixx_w_v = sum((1/12)*pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2.*(3*(0.5*(x_sect_1_w-x_sect_2_w)).^2+delta_y_w^2)+pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2.*(pos_w_y - y_w_cg).^2);  % [kg * mm^2]
    
    Iyy_w_v = sum((0.5*pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^4)+(pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2).*abs(x_sect_1_w-0.5*(x_sect_1_w-x_sect_2_w)-x_w_cg).^2); % [kg * mm^2]
    
    Izz_w_v = sum((1/12)*pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2.*(3*(0.5*(x_sect_1_w-x_sect_2_w)).^2+delta_y_w^2)+pi*delta_y_w*rho_air*(abs(x_sect_1_w-0.5*(x_sect_1_w-x_sect_2_w)-x_w_cg).^2+(pos_w_y - y_w_cg).^2));  % [kg * mm^2]
    
    %Ixy_w_v = 
    
    I_wing_v = [ Ixx_w_v 0 0; ...
                 0 Iyy_w_v 0; ...
                 0 0 Izz_w_v ]; % [ kg * mm^2 ]
           
    I_v_wing = I_wing_v+I_wing;
           
    % Return center of gravities:
    
    cg_body = [ x_b_cg; 0; 0];
    
    cg_L = [x_w_cg; -y_w_cg; 0];
    
    cg_R = [x_w_cg; y_w_cg; 0];
    
    
    % Calculate the chord-lengths of the blade sections:
    
    y_sect = zeros(3,nr_sect);
    
    y_sect(2,:) = (0.5*(w_length/nr_sect)):(w_length/nr_sect):(w_length-0.5*(w_length/nr_sect));
    
    x_sect_1_chords = interp1(y_cont_w(1:((end+1)/2)),x_cont_w(1:((end+1)/2)),y_sect(2,:));
    
    x_sect_2_chords = interp1(y_cont_w(((end+1)/2):end),x_cont_w(((end+1)/2):end),y_sect(2,:));
    
    chords = x_sect_1_chords-x_sect_2_chords;
    
    
    
    
    
end

