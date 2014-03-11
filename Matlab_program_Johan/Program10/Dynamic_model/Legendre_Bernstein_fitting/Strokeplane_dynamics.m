function [ stroke_var ] = Strokeplane_dynamics( settings, pathDB, seq_nr )


    % Return the following variables in the strokeplane reference frame:
    
    % Vn, Vt, wx, wy, wz
    
    % An, At, w_dot_x, w_dot_y, w_dot_z
    
    %----------------------------------------------------------------------
    
    stroke_var = {};
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    t = pathDB.t(start:stop);
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    
    Omega_strk = zeros(3,nr_wb);
    
    Omega_dot_strk = zeros(3,nr_wb);
    
    Vn_strk = zeros(1,nr_wb);
    
    Vt_strk = zeros(1,nr_wb);
    
    An_strk = zeros(1,nr_wb);
    
    At_strk = zeros(1,nr_wb);
    
    
    for i = 1:nr_wb
        
        a = find(isnan(pathDB.wingbeat_time(i,:,seq_nr))==0, 1, 'last' );
        
        b1 = start-1+pathDB.wingbeat_time(i,1:(a-1),seq_nr);
        
        b2 = start-1+pathDB.wingbeat_time(i,2:a,seq_nr);
        
        
        u_body = pathDB.u_body_mean(i,1,seq_nr);
        
        v_body = pathDB.v_body_mean(i,1,seq_nr);
        
        w_body = pathDB.w_body_mean(i,1,seq_nr);
        
        wx_body = pathDB.omegax_body_mean(i,1,seq_nr);
        
        wy_body = pathDB.omegay_body_mean(i,1,seq_nr);
        
        wz_body = pathDB.omegaz_body_mean(i,1,seq_nr);
        
               
        ax_body = pathDB.ax_body_mean(i,1,seq_nr);
        
        ay_body = pathDB.ay_body_mean(i,1,seq_nr);
        
        az_body = pathDB.az_body_mean(i,1,seq_nr);
        
        w_dot_x_body = mean((pathDB.b_omega1(b2,seq_nr)-pathDB.b_omega1(b1,seq_nr))./dt);
        
        w_dot_y_body = mean((pathDB.b_omega2(b2,seq_nr)-pathDB.b_omega2(b1,seq_nr))./dt);
        
        w_dot_z_body = mean((pathDB.b_omega3(b2,seq_nr)-pathDB.b_omega3(b1,seq_nr))./dt);
        
        
        % Transfer the variables to the strokeplane
        
        beta = -(55/180)*pi;    
    
        R_beta = [cos(beta) 0 -sin(beta); ...
               0 1 0; ...
               sin(beta) 0 cos(beta)]; 
           
        V_strk = R_beta'*[u_body; v_body; w_body];
        
        Omega_strk(:,i) = R_beta*[wx_body; wy_body; wz_body];
        
        A_strk = R_beta'*[ax_body; ay_body; az_body];
        
        Omega_dot_strk(:,i) = R_beta*[w_dot_x_body; w_dot_y_body; w_dot_z_body];
        
        Vn(i) = -V_strk(3);
        
        Vt(i) = sqrt(V_strk(1)^2+V_strk(2)^2);
        
        An(i) = -A_strk(3);
        
        At(i) = sqrt(A_strk(1)^2+A_strk(2)^2);
                
        
    end

    stroke_var.Omega_strk = Omega_strk;
    
    stroke_var.Omega_dot_strk = Omega_dot_strk;
    
    stroke_var.Vn = Vn;
    
    stroke_var.Vt = Vt;
    
    stroke_var.An = An;
    
    stroke_var.At = At;

end

