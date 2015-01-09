function [ F_L, F_R, c_pres_L, c_pres_R ] = Aerodynamic_Model(wing_l,y_sect, chords, joint_L, joint_R, u_b, w_b, w_L, w_R, R_b, R_L, R_R,up,down)

    % Input: body and wing kinematics.
    % Output: 2 Vectors with aerodynamic forces at the aerodynamic center
    % and the location of the center of pressure.
    
    
    % Calculate density:
    
    pres = 101.325e3; % Pressure in Pascal
    
    Temp = 297.039; % Average temperature in Kelvin
    
    R = 287.058; % Gas constant for air
    
    Rho = 1e-9*pres/(R*Temp); % Density in kg/mm^3
    
    
    nr_sect = length(y_sect);
    
    
    % Sections on left and right wing:
    
    y_sect_L = -y_sect;
    
    y_sect_R = y_sect;
    
    
    
    % Calculate the velocity at the different sections of the wings:
    
    u_L = zeros(3,nr_sect);
    
    u_R = zeros(3,nr_sect);
    
    for i = 1:nr_sect
    
        u_L(:,i) = (R_L'*(u_b + cross(w_b,joint_L + R_L*y_sect_L(:,i)))+cross(w_L,y_sect_L(:,i)));
    
        u_R(:,i) = (R_R'*(u_b + cross(w_b,joint_R + R_R*y_sect_R(:,i)))+cross(w_R,y_sect_R(:,i)));
    
    end
    
    
    
    
    % Calculate angle of attack at the sections:
    
    alfa_L = zeros(nr_sect,1);
    
    alfa_R = zeros(nr_sect,1);
    
    for i = 1:nr_sect
        
        alfa_L(i) = real(atan2(u_L(3,i),u_L(1,i)));
        
        alfa_R(i) = real(atan2(u_R(3,i),u_R(1,i)));
        
    end
    
    
    % Calculate the local lift and drag forces:
    
    deltaR = y_sect_R(2,2)-y_sect_R(2,1);
     
    Lift_L = zeros(nr_sect,1);
    
    Drag_L = zeros(nr_sect,1);
    
    Lift_R = zeros(nr_sect,1);
    
    Drag_R = zeros(nr_sect,1);
    
    F_loc_L = zeros(3,nr_sect);
    
    M_loc_L = zeros(3,nr_sect);
    
    F_loc_R = zeros(3,nr_sect);
    
    M_loc_R  =  zeros(3,nr_sect);

        
    for i = 1:nr_sect
            
        if up == 1
            
            alfa = -radtodeg(alfa_L(i));
                
            Cl = 0.225+1.58*sind(2.13*alfa-7.2);
            
            Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
            
            Lift_L(i) = 0.5*Rho*chords(i)*(u_L(1,i)^2+u_L(3,i)^2)*Cl*deltaR; % [ kg/mm^3 * mm * mm^2/s^2 * mm ]
            
            Drag_L(i) = 0.5*Rho*chords(i)*(u_L(1,i)^2+u_L(3,i)^2)*Cd*deltaR;
            
            F_loc_L(:,i) = [(Lift_L(i)*sind(alfa)-Drag_L(i)*cosd(alfa)); 0; (Lift_L(i)*cosd(alfa)+Drag_L(i)*sind(alfa))];
            
            M_loc_L(:,i) = cross(y_sect_L(:,i),F_loc_L(:,i));
            
        elseif down == 1
                
            alfa = radtodeg(alfa_L(i));
                
            Cl = 0.225+1.58*sind(2.13*alfa-7.2);
            
            Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
            
            Lift_L(i) = 0.5*Rho*chords(i)*(u_L(1,i)^2+u_L(3,i)^2)*Cl*deltaR;
            
            Drag_L(i) = 0.5*Rho*chords(i)*(u_L(1,i)^2+u_L(3,i)^2)*Cd*deltaR;
            
            F_loc_L(:,i) = [(Lift_L(i)*sind(alfa)-Drag_L(i)*cosd(alfa)); 0; -(Lift_L(i)*cosd(alfa)+Drag_L(i)*sind(alfa))];
            
            M_loc_L(:,i) = cross(y_sect_L(:,i),F_loc_L(:,i));
            
        end
            
        clear Cl Cd alfa
            
    end
    
    for i = 1:nr_sect
            
        if up == 1
            
            alfa = -radtodeg(alfa_R(i));
                
            Cl = 0.225+1.58*sind(2.13*alfa-7.2);
            
            Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
            
            Lift_R(i) = 0.5*Rho*chords(i)*(u_R(1,i)^2+u_R(3,i)^2)*Cl*deltaR;
            
            Drag_R(i) = 0.5*Rho*chords(i)*(u_R(1,i)^2+u_R(3,i)^2)*Cd*deltaR;
            
            F_loc_R(:,i) = [(Lift_R(i)*sind(alfa)-Drag_R(i)*cosd(alfa)); 0; (Lift_R(i)*cosd(alfa)+Drag_R(i)*sind(alfa))];
            
            M_loc_R(:,i) = cross(y_sect_R(:,i),F_loc_R(:,i));
            
        elseif down == 1
                
            alfa = radtodeg(alfa_R(i));
                
            Cl = 0.225+1.58*sind(2.13*alfa-7.2);
            
            Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
            
            Lift_R(i) = 0.5*Rho*chords(i)*(u_R(1,i)^2+u_R(3,i)^2)*Cl*deltaR;
            
            Drag_R(i) = 0.5*Rho*chords(i)*(u_R(1,i)^2+u_R(3,i)^2)*Cd*deltaR;
            
            F_loc_R(:,i) = [(Lift_R(i)*sind(alfa)-Drag_R(i)*cosd(alfa)); 0; -(Lift_R(i)*cosd(alfa)+Drag_R(i)*sind(alfa))];
            
            M_loc_R(:,i) = cross(y_sect_R(:,i),F_loc_R(:,i));
            
        end
            
        clear Cl Cd alfa
            
    end
%     
%     F_loc_L(1,:)
    
    % Calculate the location of the center op pressure:
    
%     sum(M_loc_L(1,:))
%     
%     sum(M_loc_L(3,:))
%     
%     sum(F_loc_L(1,:))
%     
%     sum(F_loc_R(1,:))
    
%     c_pres_L = [ 0; 0.5*(sum(M_loc_L(1,:))/sum(F_loc_L(3,:))+sum(M_loc_L(3,:))/sum(F_loc_L(1,:))); 0];
%     
%     c_pres_R = [ 0; 0.5*(sum(M_loc_R(1,:))/sum(F_loc_R(3,:))+sum(M_loc_R(3,:))/sum(F_loc_R(1,:))); 0];
 

    c_pres_L = [ 0; -0.75*wing_l; 0];
    
    c_pres_R = [ 0; 0.75*wing_l; 0];
    
    % Calculate the combined lift and drag forces at the center of
    % pressure:
    
    F_L = 1e3*sum(F_loc_L')'; % Force in Newton
    
    F_R = 1e3*sum(F_loc_R')';% Force in Newton
    

end

