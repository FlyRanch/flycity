function [ experiment ] = Aerodynamic_damping_maneuvers( nr_man, maneuvers, man_name, body_kin )


    experiment = {};

    for i = 1:nr_man
        
       temp_man = maneuvers.(char(man_name(i)));
       
       seq_nr = temp_man.seq_nr;
        
       temp_theta_LR1_avg = temp_man.a_avg.theta_LR1;
       temp_theta_LR2_avg = temp_man.a_avg.theta_LR2;
       temp_eta_LR1_avg = temp_man.a_avg.eta_LR1;
       temp_eta_LR2_avg = temp_man.a_avg.eta_LR2;
       temp_phi_LR1_avg = temp_man.a_avg.phi_LR1;
       temp_phi_LR2_avg = temp_man.a_avg.phi_LR2;
       
       temp_theta_L1_dev = temp_man.a_dev.theta_L1;
       temp_theta_L2_dev = temp_man.a_dev.theta_L2;
       temp_eta_L1_dev = temp_man.a_dev.eta_L1;
       temp_eta_L2_dev = temp_man.a_dev.eta_L2;
       temp_phi_L1_dev = temp_man.a_dev.phi_L1;
       temp_phi_L2_dev = temp_man.a_dev.phi_L2;
        
       temp_theta_R1_dev = temp_man.a_dev.theta_R1;
       temp_theta_R2_dev = temp_man.a_dev.theta_R2;
       temp_eta_R1_dev = temp_man.a_dev.eta_R1;
       temp_eta_R2_dev = temp_man.a_dev.eta_R2;
       temp_phi_R1_dev = temp_man.a_dev.phi_R1;
       temp_phi_R2_dev = temp_man.a_dev.phi_R2;
       
       f_vect = temp_man.f;
       
       down_up = temp_man.down_up;
       
       nr_wb_man = length(f_vect);
       
       for j = 1:nr_wb_man
           
           a_theta_L = [ temp_theta_LR1_avg+temp_theta_L1_dev(:,j); temp_theta_LR2_avg+temp_theta_L2_dev(:,j)];
           a_eta_L = [ temp_eta_LR1_avg+temp_eta_L1_dev(:,j); temp_eta_LR2_avg+temp_eta_L2_dev(:,j)];
           a_phi_L = [ temp_phi_LR1_avg+temp_phi_L1_dev(:,j); temp_phi_LR2_avg+temp_phi_L2_dev(:,j)];

           a_theta_R = [ temp_theta_LR1_avg+temp_theta_R1_dev(:,j); temp_theta_LR2_avg+temp_theta_R2_dev(:,j)];
           a_eta_R = [ temp_eta_LR1_avg+temp_eta_R1_dev(:,j); temp_eta_LR2_avg+temp_eta_R2_dev(:,j)];
           a_phi_R = [ temp_phi_LR1_avg+temp_phi_R1_dev(:,j); temp_phi_LR2_avg+temp_phi_R2_dev(:,j)];
           
           a_sim.theta_L(:,j) = a_theta_L;
           a_sim.eta_L(:,j) = a_eta_L;
           a_sim.phi_L(:,j) = a_phi_L;
           
           a_sim.theta_R(:,j) = a_theta_R;
           a_sim.eta_R(:,j) = a_eta_R;
           a_sim.phi_R(:,j) = a_phi_R;
           
       
       end
       
       t_vect = zeros(nr_wb_man+1,1);
          
       nr_points_vect = zeros(nr_wb_man,1);
    
       for k = 1:nr_wb_man

            nr_of_points = round((1/f_vect(k))/dt2)+1;
            nr_points_vect(k) = nr_of_points;
            t_vect(k+1) = t_vect(k)+(nr_of_points-1)*dt2;
            % t_vect(k+1) = k;
       end
            
       n_tot = sum(nr_points_vect)-nr_wb_man+1;
       
       
       nr_of_exp = length(body_kin.phi_b);
       
       
       F_exp = zeros(18,nr_of_exp,nr_wb_man);
       
       M_exp = zeros(15,nr_of_exp,nr_wb_man);
       
       
       for j = 1:nr_of_exp
            
       phi_b = body_kin.phi_b(j);
       theta_b = body_kin.theta_b(j);
       
       q_b_c  = [ sin(phi_b/2)*cos(theta_b/2); ...
                  cos(phi_b/2)*sin(theta_b/2); ...
                 -sin(phi_b/2)*sin(theta_b/2); ...
                  cos(phi_b/2)*cos(theta_b/2)];
              
       v_st_c = body_kin.v_st_c(:,j);
       a_st_c = body_kin.a_st_c(:,j);
       w_st_c = body_kin.w_st_c(:,j);
       w_dot_st_c = body_kin.w_dot_st_c(:,j);
           
       q_b = zeros(4,n_tot);
       v_st = zeros(3,n_tot);
       a_st = zeros(3,n_tot);
       w_st = zeros(3,n_tot);
       w_dot_st = zeros(3,n_tot);
           
       for k = 1:n_tot
               
         q_b(:,k) = q_b_c;
         v_st(:,k) = v_st_c;
         a_st(:,k) = a_st_c;
         w_st(:,k) = w_st_c;
         w_dot_st(:,k) = w_dot_st_c;
               
       end

       Req_f = Required_forces( settings, pathDB, a_sim, nr_points_vect, t_vect, down_up, f_vect, dt2, seq_nr, q_b, v_st, a_st, w_st, w_dot_st );
       
       Req_f
           
       t_man = Req_f.time;
       
       F_I_acc  = Req_f.F_I_acc;
       F_I_vel  = Req_f.F_I_vel;
       M_I_acc  = Req_f.M_I_acc;
       M_I_vel  = Req_f.M_I_vel;
       F_I      = Req_f.F_I;
       M_I      = Req_f.M_I;
       c_pres_L = Req_f.c_pres_L;
       c_pres_R = Req_f.c_pres_R;
       F_a_wL   = Req_f.F_a_wL;
       F_a_wR   = Req_f.F_a_wR;
       F_a_L    = Req_f.F_a_L;
       F_a_R    = Req_f.F_a_R;
       M_a_L    = Req_f.M_a_L;
       M_a_R    = Req_f.M_a_R;
       F_g      = Req_f.F_g;
       
       
        F_I_mean = Req_f.F_I_mean;
        F_I_acc_mean = Req_f.F_I_acc_mean;
        F_I_vel_mean = Req_f.F_I_vel_mean;
        F_a_L_mean = Req_f.F_a_L_mean;
        F_a_R_mean = Req_f.F_a_R_mean;
        F_g_mean = Req_f.F_g_mean;
        
        M_I_mean = Req_f.M_I_mean;
        M_I_acc_mean = Req_f.M_I_acc_mean;
        M_I_vel_mean = Req_f.M_I_vel_mean;
        M_a_L_mean = Req_f.M_a_L_mean;
        M_a_R_mean = Req_f.M_a_R_mean;

        
        for k = 1:nr_wb_man
            
            F_exp(1:3,j,k) = F_I_mean(:,k);
            F_exp(4:6,j,k) = F_I_acc_mean(:,k);
            F_exp(7:9,j,k) = F_I_vel_mean(:,k);
            F_exp(10:12,j,k) = F_a_L_mean(:,k);
            F_exp(13:15,j,k) = F_a_R_mean(:,k);
            F_exp(16:18,j,k) = F_g_mean(:,k);
            
            M_exp(1:3,j,k) = M_I_mean(:,k);
            M_exp(4:6,j,k) = M_I_acc_mean(:,k);
            M_exp(7:9,j,k) = M_I_vel_mean(:,k);
            M_exp(10:12,j,k) = M_a_L_mean(:,k);
            M_exp(13:15,j,k) = M_a_R_mean(:,k);
            
        end
        
        
        
       
%        arm_L = [0; -1; 0];
%        arm_R = [0; 1; 0];
%        
%        arm_L1 = [0.025; -1; 0];
%        arm_R1 = [0.025; 1; 0];
%        
%        arm_L2 = [-0.075; -1; 0];
%        arm_R2 = [-0.075; 1; 0];
%        
%        figure()
%        hold on
%        for k = 1:nr_points_vect(1)
%            Lwt = Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*arm_L;
%            Rwt = Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*arm_R;
%            Lwt1 = Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*arm_L1;
%            Rwt1 = Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*arm_R1;
%            Lwt2 = Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*arm_L2;
%            Rwt2 = Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*arm_R2;
%            scale = 10;
%            F_L = scale*Req_f.R_b(:,:,k)'*Req_f.R_L(:,:,k)'*Req_f.F_a_wL(:,k);
%            F_R = scale*Req_f.R_b(:,:,k)'*Req_f.R_R(:,:,k)'*Req_f.F_a_wR(:,k);
%            plot3([Lwt1(1) Lwt2(1)],[Lwt1(2) Lwt2(2)],[Lwt1(3) Lwt2(3)],'k')
%            plot3([Rwt1(1) Rwt2(1)],[Rwt1(2) Rwt2(2)],[Rwt1(3) Rwt2(3)],'k')
%            plot3(cos(0:0.05:(2*pi)),sin(0:0.05:(2*pi)),zeros(length(0:0.05:(2*pi)),1),':','Color','b')
%            plot3([-1 1],[0 0],[0 0],':','Color','b')
%            quiver3(Lwt(1),Lwt(2),Lwt(3),F_L(1),F_L(2),F_L(3),'r')
%            quiver3(Rwt(1),Rwt(2),Rwt(3),F_R(1),F_R(2),F_R(3),'r')
%        end
%        hold off
%        axis equal
%        xlabel('x')
%        ylabel('y')
%        zlabel('z')
% 
%        
%        figure()
%        hold on
%        for k = 1:nr_wb_man
%            if k == 1
%                b = 1:nr_points_vect(1);
%            else           
%                b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
%            end 
%            
%        plot(t_man(b),M_I(1,b),t_man(b),M_a_L(1,b),t_man(b),M_a_R(1,b),t_man(b),M_a_L(1,b)+M_a_R(1,b)+M_I(1,b))
%        plot(t_man(b),ones(length(b),1)*mean(M_I(1,b)),t_man(b),ones(length(b),1)*mean(M_a_L(1,b)+M_a_R(1,b)),t_man(b),ones(length(b),1)*mean(M_a_L(1,b)+M_a_R(1,b)+M_I(1,b)))
%        end
%        hold off
%        
%        figure()
%        hold on
%        for k = 1:nr_wb_man
%            if k == 1
%                b = 1:nr_points_vect(1);
%            else           
%                b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
%            end 
%            
%        plot(t_man(b),M_I(2,b),t_man(b),M_a_L(2,b),t_man(b),M_a_R(2,b),t_man(b),M_a_L(2,b)+M_a_R(2,b)+M_I(2,b))
%        plot(t_man(b),ones(length(b),1)*mean(M_I(2,b)),t_man(b),ones(length(b),1)*mean(M_a_L(2,b)+M_a_R(2,b)),t_man(b),ones(length(b),1)*mean(M_a_L(2,b)+M_a_R(2,b)+M_I(2,b)))
%        end
%        hold off
%        
%        figure()
%        hold on
%        for k = 1:nr_wb_man
%            if k == 1
%                b = 1:nr_points_vect(1);
%            else           
%                b = (sum(nr_points_vect(1:(k-1)))-(k-1)+1):(sum(nr_points_vect(1:k))-(k-1));
%            end 
%            
%        plot(t_man(b),M_I(3,b),t_man(b),M_a_L(3,b),t_man(b),M_a_R(3,b),t_man(b),M_a_L(3,b)+M_a_R(3,b)+M_I(3,b))
%        plot(t_man(b),ones(length(b),1)*mean(M_I(3,b)),t_man(b),ones(length(b),1)*mean(M_a_L(3,b)+M_a_R(3,b)),t_man(b),ones(length(b),1)*mean(M_a_L(3,b)+M_a_R(3,b)+M_I(3,b)))
%        end
%        hold off
%        
%        pause
       
       
       end
       
       experiment.(char(man_name(i))).F_exp = F_exp;
       experiment.(char(man_name(i))).M_exp = M_exp;
       
       F_exp
       
       M_exp
       
       
       pause
       
    end
    
    
    

end

