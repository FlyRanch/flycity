function [Fx,Fy,Fz,Fn,Ft,Mn,Mt,R_F,r_body,Mx,My,Mz,q_avg] = Quasi_steady_force(q1,q2,q3,q4,u_wing,v_wing,w_wing,alfa_wing,chords,deltaR,y_sect,up,down)


    % Program that calculates the quasi steady force for a given sequence
    % of wing orientations and given local angles of attack and local
    % velocities
    
    % Density is related to temperature, the average temperature is 75 deg
    % Fahrenheit (23.9 deg Celsius) with a minimum and maximum of 70 deg
    % Fahrenheit and 80 deg Fahrenheit respectively (21.11 deg Celsius and
    % 26.67 deg Celsius). The density will be derived with the perfect gas
    % law, assuming the pressure is at standard atmosphere.

    pres = 101.325e3; % Pressure in Pascal
    
    Temp = 297.039; % Average temperature in Kelvin
    
    R = 287.058; % Gas constant for air
    
    Rho = pres/(R*Temp); % Density in kg/m^3
    
    % Transfer the velocity [mm/s], chords [mm], deltaR [mm] and y_sect[mm]
    % to meter in order to generate forces in Newton:
    
    u_wing = u_wing./1000;
    v_wing = v_wing./1000;
    w_wing = w_wing./1000;
    
    chords = chords./1000;
    
    y_sect = y_sect./1000;
    
    deltaR = deltaR/1000;
    
    
    % Calculate the lift and the drag at each section for each time:
    
    N = size(u_wing,1);
    
    M = length(chords);
    
    Lift_loc = zeros(N,M);
    
    Drag_loc = zeros(N,M);
    
    F_loc = zeros(3,M,N);
    
    M_loc = zeros(3,M,N);
    

    
    for k = 1:N
        
        for l = 1:M
            
            if up == 1
            
            alfa = -radtodeg(alfa_wing(k,l));
                
            Cl = 0.225+1.58*sind(2.13*alfa-7.2);
            
            Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
            
            Lift_loc(k,l) = 0.5*Rho*chords(l)*(u_wing(k,l)^2+w_wing(k,l)^2)*Cl*deltaR;
            
            Drag_loc(k,l) = 0.5*Rho*chords(l)*(u_wing(k,l)^2+w_wing(k,l)^2)*Cd*deltaR;
            
            F_loc(:,l,k) = [(Lift_loc(k,l)*sind(alfa)-Drag_loc(k,l)*cosd(alfa)); 0; (Lift_loc(k,l)*cosd(alfa)+Drag_loc(k,l)*sind(alfa))];
            
            M_loc(:,l,k) = cross([0; y_sect(l); 0],F_loc(:,l,k));
            
            elseif down == 1
                
            alfa = radtodeg(alfa_wing(k,l));
                
            Cl = 0.225+1.58*sind(2.13*alfa-7.2);
            
            Cd = 1.92-1.55*cosd(2.04*alfa-9.82);
            
            Lift_loc(k,l) = 0.5*Rho*chords(l)*(u_wing(k,l)^2+w_wing(k,l)^2)*Cl*deltaR;
            
            Drag_loc(k,l) = 0.5*Rho*chords(l)*(u_wing(k,l)^2+w_wing(k,l)^2)*Cd*deltaR;
            
            F_loc(:,l,k) = [(Lift_loc(k,l)*sind(alfa)-Drag_loc(k,l)*cosd(alfa)); 0; -(Lift_loc(k,l)*cosd(alfa)+Drag_loc(k,l)*sind(alfa))];
            
            M_loc(:,l,k) = cross([0; y_sect(l); 0],F_loc(:,l,k));
            
            end
            
            clear Cl Cd alfa
            
        end
        
    end



    
    % Transfer all the lift and drag forces to the body frame of reference
    % and average the local lift and drag forces at the center of pressure:
    
    Fn = nan(1,50); %Sum of the spanwise normal force.
    
    Ft = nan(1,50); %Sum of the spanwise tangential force.
    
    Mn = nan(1,50); %Sum of the spanwise normal moment.
    
    Mt = nan(1,50); %Sum of the spanwize tangential force.
    
    R_F = nan(1,50);
    
    for k = 1:N
               
        Fn(k) = sum(F_loc(3,:,k));
        
        Ft(k) = sum(F_loc(1,:,k));
        
        Mn(k) = sum(M_loc(1,:,k));
        
        Mt(k) = sum(M_loc(3,:,k));
        
        R_F(k) = sign(y_sect(1))*1000*sqrt(Mn(k)^2+Mt(k)^2)/sqrt(Fn(k)^2+Ft(k)^2);
        
    end
    
    % Transfer the local forces and moments to the body frame of reference
    % at the joint location:
    
    F_joint = zeros(3,N);
    
    M_joint = zeros(3,N);
    
    R_joint = zeros(3,N);
    
    for k = 1:N
        
        F_loc_body = zeros(3,M);
        
        %M_loc_body = zeros(3,M);
        
        DCM = quat2matNEW([q1(k) q2(k) q3(k) q4(k)]);
        
        for l = 1:M
                        
            F_loc_body(:,l) = DCM*F_loc(:,l,k);
            
            %M_loc_body(:,l) = DCM'*M_loc(:,l,k);
            
                       
        end
        
        R_joint(:,k) = DCM*[0; R_F(k); 0];
        
        F_joint(1,k) = sum(F_loc_body(1,:));
        F_joint(2,k) = sum(F_loc_body(2,:));
        F_joint(3,k) = sum(F_loc_body(3,:));
        
%         F_joint(:,k) = sum(F_loc_body,2);
        
        %M_joint(:,k) = sum(M_loc_body,2);
        
        M_joint(:,k) = DCM*[sum(M_loc(1,:,k)); sum(M_loc(2,:,k)); sum(M_loc(3,:,k))];
%         M_joint(2,k) = sum(M_loc(2,:,k));
%         M_joint(3,k) = sum(M_loc(3,:,k));
%         
        clear F_loc_body M_loc_body DCM
        
    end
    
%     figure()
%     plot(F_joint(1,1:N))
%     hold on
%     plot(F_joint(2,1:N),'r')
%     plot(F_joint(3,1:N),'g')
%     hold off
%     
%     figure()
%     plot(M_joint(1,1:N))
%     hold on
%     plot(M_joint(2,1:N),'r')
%     plot(M_joint(3,1:N),'g')
%     hold off
        
    %Calculate the length and orientation of the arm necessary to generate
    
        Fx = F_joint(1,:);
        
        Fy = F_joint(2,:);
        
        Fz = F_joint(3,:);
        
        Mx = M_joint(1,:);
         
        My = M_joint(2,:);
        
        Mz = M_joint(3,:);
        
        r_body = mean(R_joint,2);
%         
%         M_t = cross(r_body./1000,[mean(;Fy;Fz]);
%         
%         Mx = M_t(1);
%          
%         My = M_t(2);
%         
%         Mz = M_t(3);

    
    % Determine the average orientation of wing during the given time
    % sequence:
    
    % Initial estimate:
    
    m = [q1(1); q2(1); q3(1); q4(1)];
    
    err = 1;
    
    while abs(err)<0.001
        
        err_t = ones(M,1);
        
        for k = 1:M
            
            m_inv = [-m(1); -m(2); -m(3); m(4)]./(m(1)^2+m(2)^2+m(3)^2+m(4)^2);
            
            m_inv = m_inv./norm(m_inv);
            
            M_inv = [m_inv(4) m_inv(3) -m_inv(2) m_inv(1); ...
                     -m_inv(3) m_inv(4) m_inv(1) m_inv(2); ...
                     m_inv(2) -m_inv(1) m_inv(4) m_inv(3); ...
                     -m_inv(1) -m_inv(2) -m_inv(3) m_inv(4)];

            M_q = M_inv*[q1(k); q2(k); q3(k); q4(k)];
            
            M_q = M_q./norm(M_q);
                 
            err_t(k) = log(M_q);
            
            clear m_inv M_inv M_q
            
        end
        
        err = sum(err_t);
        
        m = m*exp(err);
        
        m = m./norm(m);
        
        clear err_t
        
    end
    
    q_avg = [m(1); m(2); m(3); m(4)];



end

