function Wingkinematics_Legendre(settings, pathDB)


    % Use a Legendre polynomial fit to describe the wingkinematics in a
    % minimal number of parameters:
    
    %     n_pol_theta = 13; % Order of used polynomials
%     n_pol_eta = 16; % Order of used polynomials
%     n_pol_phi = 8; % Order of used polynomials
%     
    n_pol_theta = 12; % Order of used polynomials
    n_pol_eta = 14; % Order of used polynomials
    n_pol_phi = 10; % Order of used polynomials
%     n_pol_theta = 9; % Order of used polynomials
%     n_pol_eta = 9; % Order of used polynomials
%     n_pol_phi = 7; % Order of used polynomials


    n_pol = {};
    
    n_pol.theta = n_pol_theta;
    
    n_pol.eta = n_pol_eta;
    
    n_pol.phi = n_pol_phi;
    

    
    dt = pathDB.t(2)-pathDB.t(1);
    
    nr_of_seq = size(pathDB.x,2);
    
    Maneuver_Matrix = zeros(1,13);
    
    a_sym_maneuver  = {};
    
    a_dev_maneuver = {};
    
    
    a_sym_maneuver.theta_L = zeros(1,n_pol_theta+1);
    
    a_sym_maneuver.theta_R = zeros(1,n_pol_theta+1);
    
    a_sym_maneuver.eta_L = zeros(1,n_pol_eta+1);
    
    a_sym_maneuver.eta_R = zeros(1,n_pol_eta+1);
    
    a_sym_maneuver.phi_L = zeros(1,n_pol_phi+1);
    
    a_sym_maneuver.phi_R = zeros(1,n_pol_phi+1);
    
    
    
    a_dev_maneuver.theta_L = zeros(1,n_pol_theta+1);
    
    a_dev_maneuver.theta_R = zeros(1,n_pol_theta+1);
    
    a_dev_maneuver.eta_L = zeros(1,n_pol_eta+1);
    
    a_dev_maneuver.eta_R = zeros(1,n_pol_eta+1);
    
    a_dev_maneuver.phi_L = zeros(1,n_pol_phi+1);
    
    a_dev_maneuver.phi_R = zeros(1,n_pol_phi+1);
    
    
    for i =1:nr_of_seq
    
    seq_nr = i
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );


    [a_fit,a_avg,f_avg,down_up,trigger_wb,ratio_1_2,ratio_1_2_avg] = Standard_wingbeat( settings, pathDB, i, n_pol_theta, n_pol_eta, n_pol_phi);

   
    
    if trigger_wb >= 3 && trigger_wb <(nr_wb-1)


    [a_sym, a_dev] = Maneuver_wingbeat(a_fit,a_avg,n_pol_theta,n_pol_eta,n_pol_phi,trigger_wb,ratio_1_2,ratio_1_2_avg);
    
    
    
    [ stroke_var ] = Strokeplane_dynamics( settings, pathDB, i );

    
    
    dev_theta_L = zeros(100*nr_wb,1);
    
    dev_theta_R = zeros(100*nr_wb,1);
    
    dev_eta_L = zeros(100*nr_wb,1);
    
    dev_eta_R = zeros(100*nr_wb,1);
    
    dev_phi_L = zeros(100*nr_wb,1);
    
    dev_phi_R = zeros(100*nr_wb,1);
    
    dev_theta_dot_L = zeros(100*nr_wb,1);
    
    dev_theta_dot_R = zeros(100*nr_wb,1);
    
    dev_eta_dot_L = zeros(100*nr_wb,1);
    
    dev_eta_dot_R = zeros(100*nr_wb,1);
    
    dev_phi_dot_L = zeros(100*nr_wb,1);
    
    dev_phi_dot_R = zeros(100*nr_wb,1);
    
    time_loc = pathDB.wingbeat_time(1:nr_wb,1,seq_nr);
    
    time = pathDB.t(start-1+time_loc);
    
    time2 = zeros(100*nr_wb,1);
    
  
    for k = 1:nr_wb
        
        if k < nr_wb
    
            [ ~ , X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(k), 100, time(k), time(k+1), 0 );

            [ ~ , X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2(k), 100, time(k), time(k+1), 0 );

            [ t_temp , X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2(k), 100, time(k), time(k+1), 0 );

            time2(((k-1)*100+1):(k*100)) = t_temp;

            [ ~ , X_dot_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(k), 100, time(k), time(k+1), 1 );

            [ ~ , X_dot_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2(k), 100, time(k), time(k+1), 1 );

            [ ~ , X_dot_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2(k), 100, time(k), time(k+1), 1 );

            dev_theta_L(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_R(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_L(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_R(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_L(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_R(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];

            dev_theta_dot_L(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_dot_R(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_dot_L(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_dot_R(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_dot_L(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_dot_R(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];
        
        elseif k == nr_wb
            
            t_end = find(isnan(pathDB.wingbeat_time(nr_wb,:,seq_nr))==0, 1, 'last' );
            
            [ ~ , X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(k), 100, time(k), time(k)+(t_end-1)*dt, 0 );

            [ ~ , X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2(k), 100, time(k), time(k)+(t_end-1)*dt, 0 );

            [ t_temp , X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2(k), 100, time(k), time(k)+(t_end-1)*dt, 0 );

            time2(((k-1)*100+1):(k*100)) = t_temp;

            [ ~ , X_dot_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, ratio_1_2(k), 100, time(k), time(k)+(t_end-1)*dt, 1 );

            [ ~ , X_dot_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, ratio_1_2(k), 100, time(k), time(k)+(t_end-1)*dt, 1 );

            [ ~ , X_dot_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, ratio_1_2(k), 100, time(k), time(k)+(t_end-1)*dt, 1 );

            dev_theta_L(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_R(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_L(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_R(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_L(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_R(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];

            dev_theta_dot_L(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_dot_R(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_dot_L(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_dot_R(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_dot_L(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_dot_R(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];    
            
        end

    end
    
    dev_data = {};
    
    dev_data.theta_L = dev_theta_L;
    
    dev_data.theta_R = dev_theta_R;
    
    dev_data.eta_L = dev_eta_L;
    
    dev_data.eta_R = dev_eta_R;
    
    dev_data.phi_L = dev_phi_L;
    
    dev_data.phi_R = dev_phi_R;
    
    
    c_var = 1.2; % threshold maneuvering flight ( c_var times variance steady flight ).
    
    maneuver_wb = select_maneuver( dev_data,trigger_wb, nr_wb, c_var );
    
    maneuver_wb.id
    
    
    nr_of_maneuvers = length(maneuver_wb.id);
    
    temp_maneuver_matrix = zeros(nr_of_maneuvers,13);
    
    for j = 1:nr_of_maneuvers
        
        temp_maneuver_matrix(j,:) = [ maneuver_wb.id(j) i ratio_1_2(j) stroke_var.Omega_strk(:,maneuver_wb.id(j))' stroke_var.Omega_dot_strk(:,maneuver_wb.id(j))' ...
                                      stroke_var.Vn(maneuver_wb.id(j)) stroke_var.Vt(maneuver_wb.id(j)) stroke_var.An(maneuver_wb.id(j)) stroke_var.At(maneuver_wb.id(j))];
        
    end
    
    
    
  
    
    
    if i == 1
        
        Maneuver_Matrix = temp_maneuver_matrix;
        
        a_sym_maneuver.theta_L = [a_sym.theta_L1; a_sym.theta_L2];

        a_sym_maneuver.theta_R = [a_sym.theta_R1; a_sym.theta_R2];

        a_sym_maneuver.eta_L = [a_sym.eta_L1; a_sym.eta_L2];

        a_sym_maneuver.eta_R = [a_sym.eta_R1; a_sym.eta_R2];

        a_sym_maneuver.phi_L = [a_sym.phi_L1; a_sym.phi_L2];

        a_sym_maneuver.phi_R = [a_sym.phi_R1; a_sym.phi_R2];
       
    
        a_dev_maneuver.theta_L = [a_dev.theta_L1; a_dev.theta_L2];

        a_dev_maneuver.theta_R = [a_dev.theta_R1; a_dev.theta_R2];

        a_dev_maneuver.eta_L = [a_dev.eta_L1; a_dev.eta_L2];

        a_dev_maneuver.eta_R = [a_dev.eta_R1; a_dev.eta_R2];

        a_dev_maneuver.phi_L = [a_dev.phi_L1; a_dev.phi_L2];

        a_dev_maneuver.phi_R = [a_dev.phi_R1; a_dev.phi_R2];
        
    else
        
        Maneuver_Matrix = [Maneuver_Matrix; temp_maneuver_matrix];
        
        a_sym_maneuver.theta_L = [a_sym_maneuver.theta_L [a_sym.theta_L1; a_sym.theta_L2]];

        a_sym_maneuver.theta_R = [a_sym_maneuver.theta_R [a_sym.theta_R1; a_sym.theta_R2]];

        a_sym_maneuver.eta_L = [a_sym_maneuver.eta_L [a_sym.eta_L1; a_sym.eta_L2]];

        a_sym_maneuver.eta_R = [a_sym_maneuver.eta_R [a_sym.eta_R1; a_sym.eta_R2]];

        a_sym_maneuver.phi_L = [a_sym_maneuver.phi_L [a_sym.phi_L1; a_sym.phi_L2]];

        a_sym_maneuver.phi_R = [a_sym_maneuver.phi_R [a_sym.phi_R1; a_sym.phi_R2]];
       
    
        a_dev_maneuver.theta_L = [a_dev_maneuver.theta_L [a_dev.theta_L1; a_dev.theta_L2]];

        a_dev_maneuver.theta_R = [a_dev_maneuver.theta_R [a_dev.theta_R1; a_dev.theta_R2]];

        a_dev_maneuver.eta_L = [a_dev_maneuver.eta_L [a_dev.eta_L1; a_dev.eta_L2]];

        a_dev_maneuver.eta_R = [a_dev_maneuver.eta_R [a_dev.eta_R1; a_dev.eta_R2]];

        a_dev_maneuver.phi_L = [a_dev_maneuver.phi_L [a_dev.phi_L1; a_dev.phi_L2]];

        a_dev_maneuver.phi_R = [a_dev_maneuver.phi_R [a_dev.phi_R1; a_dev.phi_R2]];
        
    end
    
    
    
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(time2,radtodeg(dev_theta_L),'r',time2,radtodeg(dev_theta_R),'g');
%     subplot(3,1,2); plot(time2,radtodeg(dev_eta_L),'r',time2,radtodeg(dev_eta_R),'g');
%     subplot(3,1,3); plot(time2,radtodeg(dev_phi_L),'r',time2,radtodeg(dev_phi_R),'g');
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(time2,radtodeg(dev_theta_L)-radtodeg(dev_theta_R));
%     subplot(3,1,2); plot(time2,radtodeg(dev_eta_L)-radtodeg(dev_eta_R));
%     subplot(3,1,3); plot(time2,radtodeg(dev_phi_L)-radtodeg(dev_phi_R));
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(time2,(radtodeg(dev_theta_L)+radtodeg(dev_theta_R))/2);
%     subplot(3,1,2); plot(time2,(radtodeg(dev_eta_L)+radtodeg(dev_eta_R))/2);
%     subplot(3,1,3); plot(time2,(radtodeg(dev_phi_L)+radtodeg(dev_phi_R))/2);
%     hold off
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(stroke_var.Omega_strk(1,:))
%     subplot(3,1,2); plot(stroke_var.Omega_strk(2,:))
%     subplot(3,1,3); plot(stroke_var.Omega_strk(3,:))
%     hold off
    

    

    end
    
    end

    
    a_sym_maneuver
    
    a_dev_maneuver


Maneuver_Matrix


% [w_x_strk_sorted, w_x_strk_id] = sort(Maneuver_Matrix(:,4));
% 
% [w_y_strk_sorted, w_y_strk_id] = sort(Maneuver_Matrix(:,5));
% 
% [w_z_strk_sorted, w_z_strk_id] = sort(Maneuver_Matrix(:,6));
% 
% [w_dot_x_strk_sorted, w_dot_x_strk_id] = sort(Maneuver_Matrix(:,7));
% 
% [w_dot_y_strk_sorted, w_dot_y_strk_id] = sort(Maneuver_Matrix(:,8));
% 
% [w_dot_z_strk_sorted, w_dot_z_strk_id] = sort(Maneuver_Matrix(:,9));
% 
% [Vn_strk_sorted, Vn_strk_id] = sort(Maneuver_Matrix(:,10));
% 
% [Vt_strk_sorted, Vt_strk_id] = sort(Maneuver_Matrix(:,11));
% 
% [An_strk_sorted, An_strk_id] = sort(Maneuver_Matrix(:,12));
% 
% [At_strk_sorted, At_strk_id] = sort(Maneuver_Matrix(:,13));
% 
% 
% nr_of_man_wb = length(Maneuver_Matrix(:,1));
% 
% theta_L_surf_sym = zeros(nr_of_man_wb,100);
% 
% theta_R_surf_sym = zeros(nr_of_man_wb,100);
% 
% theta_L_surf_dev = zeros(nr_of_man_wb,100);
% 
% theta_R_surf_dev = zeros(nr_of_man_wb,100);
% 
% 
% for k = 1:nr_of_man_wb
%    
%     [ ~ , X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, Maneuver_Matrix(k,3), 100, 0, 1, 0 );
%     
%     theta_L_surf_sym(k,:) = X_theta*a_sym_maneuver.theta_L(:,k);
%     
%     theta_R_surf_sym(k,:) = X_theta*a_sym_maneuver.theta_R(:,k);
%     
%     theta_L_surf_dev(k,:) = X_theta*a_dev_maneuver.theta_L(:,k);
%     
%     theta_R_surf_dev(k,:) = X_theta*a_dev_maneuver.theta_R(:,k);
%     
%     
% end
% 
% figure()
% surf(theta_L_surf_sym(w_x_strk_id,:))
% 
% figure()
% surf(theta_R_surf_sym(w_x_strk_id,:))
% 
% figure()
% surf(theta_L_surf_dev(w_x_strk_id,:))
% 
% figure()
% surf(theta_R_surf_dev(w_x_strk_id,:))
% 
% figure()
% surf(theta_L_surf_sym(w_y_strk_id,:))
% 
% figure()
% surf(theta_R_surf_sym(w_y_strk_id,:))
% 
% figure()
% surf(theta_L_surf_dev(w_y_strk_id,:))
% 
% figure()
% surf(theta_R_surf_dev(w_y_strk_id,:))
% 
% figure()
% surf(theta_L_surf_sym(w_z_strk_id,:))
% 
% figure()
% surf(theta_R_surf_sym(w_z_strk_id,:))
% 
% figure()
% surf(theta_L_surf_dev(w_z_strk_id,:))
% 
% figure()
% surf(theta_R_surf_dev(w_z_strk_id,:))
%     
% figure()
% surf(theta_L_surf_sym(w_dot_x_strk_id,:))
% 
% figure()
% surf(theta_R_surf_sym(w_dot_x_strk_id,:))
% 
% figure()
% surf(theta_L_surf_dev(w_dot_x_strk_id,:))
% 
% figure()
% surf(theta_R_surf_dev(w_dot_x_strk_id,:))
% 
% figure()
% surf(theta_L_surf_sym(w_dot_y_strk_id,:))
% 
% figure()
% surf(theta_R_surf_sym(w_dot_y_strk_id,:))
% 
% figure()
% surf(theta_L_surf_dev(w_dot_y_strk_id,:))
% 
% figure()
% surf(theta_R_surf_dev(w_dot_y_strk_id,:))
% 
% figure()
% surf(theta_L_surf_sym(w_dot_z_strk_id,:))
% 
% figure()
% surf(theta_R_surf_sym(w_dot_z_strk_id,:))
% 
% figure()
% surf(theta_L_surf_dev(w_dot_z_strk_id,:))
% 
% figure()
% surf(theta_R_surf_dev(w_dot_z_strk_id,:))


