function [cg_body,cg_L,cg_R,I_body,I_wing,I_v_wing,m_w,m_vw,y_sect,chords,area_w,x_LE] = comp_Inertia(rho_air,rho_cuticle,h_wing,b_length,w_length,mass_fly,body_model,wing_model,J_pos_R,nr_sect)


    % Set number of elements, position of the elements and in x and y direction for body and wing
    
    sect_nr_b_x = 100;
    sect_nr_w_x = 50;
    sect_nr_w_y = 100;
        
    x_cont_b = [body_model.x_mod(1:15,7); body_model.x_mod(17:21,7)];
    y_cont_b = [body_model.y_mod(1:15,7); body_model.y_mod(17:21,7)];
    x_cont_w = -J_pos_R(1) + wing_model.x_mod(:,1);
    y_cont_w = -J_pos_R(2) + wing_model.y_mod(:,1);
    
    
    clear temp_x_cont_b temp_y_cont_b temp_x_cont_w temp_y_cont_w
    
    b_width = 2*max(y_cont_b);
    
    w_width = max(x_cont_w)-min(x_cont_w);
    
    delta_x_b = b_length/sect_nr_b_x;
    delta_x_w = w_width/sect_nr_w_x;
    delta_y_w = w_length/sect_nr_w_y;
    
    pos_b_x = min(x_cont_b)+((0.5/(sect_nr_b_x)):(1/(sect_nr_b_x)):(1-(0.5/(sect_nr_b_x)))).*b_length;
    pos_w_x = min(x_cont_w)+((0.5/(sect_nr_w_x)):(1/(sect_nr_w_x)):(1-(0.5/(sect_nr_w_x)))).*w_width;
    pos_w_y = ((0.5/(sect_nr_w_y)):(1/(sect_nr_w_y)):(1-(0.5/(sect_nr_w_y)))).*w_length;

%     pos_w_x = [min(x_cont_w) min(x_cont_w)+((0.5/(sect_nr_w_x)):(1/(sect_nr_w_x)):(1-(0.5/(sect_nr_w_x)))).*w_width min(x_cont_w)+w_width];
%     pos_w_y = [ 0 ((0.5/(sect_nr_w_y)):(1/(sect_nr_w_y)):(1-(0.5/(sect_nr_w_y)))).*w_length w_length ];
    
    [Xw, Yw] = meshgrid(pos_w_x, pos_w_y);
    
%     figure()
%     hold on
%     plot(x_cont_b,y_cont_b)
%     plot(pos_b_x,zeros(sect_nr_b_x,1),'o')
%     axis equal
%     hold off
%     
%     figure()
%     hold on
%     plot(x_cont_w,y_cont_w)
% %     plot(pos_w_x,zeros(sect_nr_w_x,1),'o')
% %     plot(zeros(sect_nr_w_y,1),pos_w_y,'o')
%     plot(Xw,Yw,'o')
%     axis equal
%     hold off

    
    
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
    
%     figure()
%     hold on
%     plot(x_sect_1_w,pos_w_y,'o','Color','r')
%     plot(x_sect_2_w,pos_w_y,'o','Color','b')
%     plot(x_cont_w,y_cont_w,'k')
%     axis equal
%     hold off
%     
%     figure()
%     hold on
%     plot(pos_w_x,y_sect_1_w,'o','Color','r')
%     plot(pos_w_x,y_sect_2_w,'o','Color','b')
%     plot(x_cont_w,y_cont_w,'k')
%     axis equal
%     hold off
    
    Zw = zeros(sect_nr_w_y,sect_nr_w_x);
    Aw = delta_x_w*delta_y_w;
    
    
    for j = 1:sect_nr_w_y
        for i = 1:sect_nr_w_x
            if x_sect_1_w(j) >= Xw(j,i) && Xw(j,i) >= x_sect_2_w(j)
                if y_sect_1_w(i) >= Yw(j,i) && Yw(j,i) >= y_sect_2_w(i)
                    Zw(j,i) = 1;
                end
            end
        end
    end

%     size(Xw)
%     size(Yw)
%     
%     Zw = zeros(sect_nr_w_y+2,sect_nr_w_x+2);
%     Aw = delta_x_w*delta_y_w;
%     
%     
%     for j = 2:(sect_nr_w_y+1)
%         for i = 2:(sect_nr_w_x+1)
%             if x_sect_1_w(j) >= Xw(j,i) && Xw(j,i) >= x_sect_2_w(j)
%                 if y_sect_1_w(i) >= Yw(j,i) && Yw(j,i) >= y_sect_2_w(i)
%                     Zw(j,i) = 1;
%                 end
%             end
%         end
%     end
%     
%     size(Zw)
% 
%     figure()
%     hold on
%     surf(Xw,Yw,Zw)
%     axis equal
%     hold off
    
    
    
    % Determine volume body:
    
    vol_b = sum(2*pi*0.5*y_sect_b.*y_sect_b*delta_x_b); % [mm^3]
    
    
    % Determine area wing:
    
%     area_w = sum((x_sect_1_w-x_sect_2_w)*delta_y_w); % [mm^2]
    
    area_w = sum(sum(Aw.*Zw));
    
    
    % Determine virtual air volume wing:
    
    vol_v_w = sum(pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w);    
    
    % Calculate wing and body mass:
       
    m_w = w_length*area_w*rho_cuticle*h_wing;   % [kg]
    
    m_b = mass_fly - 2*m_w;   % [kg]
    
    m_v_w = vol_v_w*rho_air;   % [kg]
    
    m_vw = m_w+m_v_w;          % [kg]
    
   
    % Determine body density:
    
    rho_body = m_b/vol_b; % density in [kg/mm^3]
    

    % Determine wing density (2D):
    
    rho_w = m_w/area_w; % density (2D) in [kg/mm^2]
    
    
    % Determine the center of gravity of the body:
    
    x_b_cg = sum((y_sect_b.*delta_x_b).*pos_b_x)/sum(y_sect_b.*delta_x_b);

%     x_b_cg = 0; % Set to zero as the body model misses the abodomen often.
    
    % Determine the center of gravity of the wing:
    
    x_w_cg = sum(((y_sect_1_w-y_sect_2_w)*delta_x_w).*pos_w_x)/area_w;
    y_w_cg = sum(((x_sect_1_w-x_sect_2_w)*delta_y_w).*pos_w_y)/area_w;
    
    % Determine the center of gravity of the wing
    
    x_vw_cg = sum((pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w).*((x_sect_1_w+x_sect_2_w)/2))/vol_v_w;
    y_vw_cg = sum((pi*(0.5*(x_sect_1_w-x_sect_2_w)).^2*delta_y_w).*pos_w_y)/vol_v_w;
    
    % Compute combined center of gravity:
    
    x_cg_wing = (x_w_cg*m_w+x_vw_cg*m_v_w)/m_vw;
    y_cg_wing = (y_w_cg*m_w+y_vw_cg*m_v_w)/m_vw;
    
    % Determine the inertia tensor of the body:
    
    Ixx_b = sum(0.5*pi*rho_body*delta_x_b*y_sect_b.^4); % [kg * mm^2]
    
    Iyy_b = sum((1/12)*pi*rho_body*delta_x_b*y_sect_b.^2.*(3*y_sect_b.^2+delta_x_b^2)+pi*rho_body*delta_x_b*y_sect_b.^2.*abs(pos_b_x-x_b_cg).^2); % [kg * mm^2]
    
    Izz_b = Iyy_b; % [kg * mm^2]
    
        
    I_body = [ Ixx_b 0 0; ...
               0 Iyy_b 0; ...
               0 0 Izz_b];
           
    
%     % Determine the inertia tensor of the wing:
%     
%     Ixx_w = sum(delta_y_w*rho_w*(x_sect_1_w-x_sect_2_w).*(pos_w_y - y_w_cg).^2); % [kg * mm^2]
%     Iyy_w = sum(delta_x_w*rho_w*(y_sect_1_w-y_sect_2_w).*(pos_w_x - x_w_cg).^2); % [kg * mm^2]
%     Izz_w = sum((1/12)*delta_y_w*rho_w*(x_sect_1_w-x_sect_2_w).*((x_sect_1_w-x_sect_2_w).^2+delta_y_w.^2)+delta_y_w*rho_w*(x_sect_1_w-x_sect_2_w).*((pos_w_y - y_w_cg).^2+(pos_w_x - x_w_cg).^2)); % [kg * mm^2]
% 
%     I_wing = [ Ixx_w 0 0; ...
%                0 Iyy_w 0; ...
%                0 0 Izz_w ]; % [ kg * mm^2 ]

    % Determine the inertia tensor of the wing:
    
    Ixx_w_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    Iyy_w_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    Izz_w_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    Ixy_w_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    
    for j = 1:sect_nr_w_y
        for i = 1:sect_nr_w_x
            Ixx_w_temp(j,i) = rho_w*Aw*Zw(j,i)*((pos_w_y(j) - y_cg_wing)^2);
            Iyy_w_temp(j,i) = rho_w*Aw*Zw(j,i)*((pos_w_x(i) - x_cg_wing)^2);
            Izz_w_temp(j,i) = rho_w*Aw*Zw(j,i)*((pos_w_y(j) - y_cg_wing)^2+(pos_w_x(i) - x_cg_wing)^2);
            Ixy_w_temp(j,i) = rho_w*Aw*Zw(j,i)*(pos_w_x(i) - x_cg_wing)*(pos_w_y(j) - y_cg_wing);
        end
    end
    
    Ixx_w = sum(sum(Ixx_w_temp));
    Iyy_w = sum(sum(Iyy_w_temp));
    Izz_w = sum(sum(Izz_w_temp));
    Ixy_w = sum(sum(Ixy_w_temp));
    
    I_wing = [ Ixx_w -Ixy_w 0; ...
               -Ixy_w Iyy_w 0; ...
               0 0 Izz_w ]; % [ kg * mm^2 ]
           
           
%     % Determine the inertia tensor of the wing including virtual wing mass:
%     
%     Ixx_w_v = sum((1/12)*pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2.*(3*(0.5*(x_sect_1_w-x_sect_2_w)).^2+delta_y_w^2)+pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2.*(pos_w_y - y_w_cg).^2);  % [kg * mm^2]
%     Iyy_w_v = sum((0.5*pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^4)+(pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2).*abs(x_sect_1_w-0.5*(x_sect_1_w-x_sect_2_w)-x_w_cg).^2); % [kg * mm^2]
%     Izz_w_v = sum((1/12)*pi*delta_y_w*rho_air*(0.5*(x_sect_1_w-x_sect_2_w)).^2.*(3*(0.5*(x_sect_1_w-x_sect_2_w)).^2+delta_y_w^2)+pi*delta_y_w*rho_air*(abs(x_sect_1_w-0.5*(x_sect_1_w-x_sect_2_w)-x_w_cg).^2+(pos_w_y - y_w_cg).^2));  % [kg * mm^2]
%     
%     I_wing_v = [ Ixx_w_v 0 0; ...
%                  0 Iyy_w_v 0; ...
%                  0 0 Izz_w_v ]; % [ kg * mm^2 ]

    Zh = zeros(sect_nr_w_y,sect_nr_w_x);
    
    
    for j = 1:sect_nr_w_y
        for i = 1:sect_nr_w_x
            if Zw(j,i) == 1
                
                Zh(j,i) = sqrt(((x_sect_1_w(j)-x_sect_2_w(j))/2)^2-(pos_w_x(i)-(x_sect_1_w(j)-((x_sect_1_w(j)-x_sect_2_w(j))/2)))^2);
                
            end
        end
    end
    
% 
%     Zh = zeros(sect_nr_w_y+2,sect_nr_w_x+2);
%     
%     
%     for j = 1:(sect_nr_w_y+2)
%         for i = 1:(sect_nr_w_x+2)
%             if Zw(j,i) == 1
%                 
%                 Zh(j,i) = sqrt(((x_sect_1_w(j)-x_sect_2_w(j))/2)^2-(pos_w_x(i)-(x_sect_1_w(j)-((x_sect_1_w(j)-x_sect_2_w(j))/2)))^2);
%                 
%             end
%         end
%     end
%     
%     
%     figure()
%     hold on
%     surf(Xw,Yw,Zh)
%     surf(Xw,Yw,-Zh)
%     axis equal
%     hold off
%     
%     pause
    
    Ixx_vw_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    Iyy_vw_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    Izz_vw_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    Ixy_vw_temp = zeros(sect_nr_w_y,sect_nr_w_x);
    
    for j = 1:sect_nr_w_y
        for i = 1:sect_nr_w_x
            Ixx_vw_temp(j,i) = rho_air*2*Aw*Zh(j,i)*((pos_w_y(j) - y_cg_wing)^2)+(1/12)*rho_air*2*Aw*Zh(j,i)*((2*Zh(j,i))^2+delta_y_w^2);
            Iyy_vw_temp(j,i) = rho_air*2*Aw*Zh(j,i)*((pos_w_x(i) - x_cg_wing)^2)+(1/12)*rho_air*2*Aw*Zh(j,i)*((2*Zh(j,i))^2+delta_x_w^2);
            Izz_vw_temp(j,i) = rho_air*2*Aw*Zh(j,i)*((pos_w_y(j) - y_cg_wing)^2+(pos_w_x(i) - x_cg_wing)^2);
            Ixy_vw_temp(j,i) = rho_air*2*Aw*Zh(j,i)*(pos_w_x(i) - x_cg_wing)*(pos_w_y(j) - y_cg_wing);
        end
    end
    
    Ixx_vw = sum(sum(Ixx_vw_temp));
    Iyy_vw = sum(sum(Iyy_vw_temp));
    Izz_vw = sum(sum(Izz_vw_temp));
    Ixy_vw = sum(sum(Ixy_vw_temp));
    
    I_v_wing = [ Ixx_vw -Ixy_vw 0; ...
                -Ixy_vw Iyy_vw 0; ...
                 0 0 Izz_vw ]; % [ kg * mm^2 ]
                        
    % Return center of gravities:
    
    cg_body = [ x_b_cg; 0; 0];
    cg_L = [x_cg_wing; -y_cg_wing; 0];
    cg_R = [x_cg_wing; y_cg_wing; 0];

    I_v_wing = I_v_wing+I_wing;
    
    
    % Calculate the chord-lengths of the blade sections:
    
    y_sect = zeros(3,nr_sect);
    y_sect(2,:) = (0.5*(w_length/nr_sect)):(w_length/nr_sect):(w_length-0.5*(w_length/nr_sect));
    x_sect_1_chords = interp1(y_cont_w(1:((end+1)/2)),x_cont_w(1:((end+1)/2)),y_sect(2,:));
    x_sect_2_chords = interp1(y_cont_w(((end+1)/2):end),x_cont_w(((end+1)/2):end),y_sect(2,:));

    chords = x_sect_1_chords-x_sect_2_chords;
    
    x_LE = x_sect_1_chords;


end

