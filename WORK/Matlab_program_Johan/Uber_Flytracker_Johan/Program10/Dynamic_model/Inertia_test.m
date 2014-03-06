function Inertia_test( wing_kin, Body, time )

    % Plot the inertial forces and moments for a given set of wing
    % kinematics:
    
    
    % Sequence length:
    
    n = length(time);
    

    
    % Execute different test-case
    
    test_case_nr = 3;
    
    switch test_case_nr
        
        case 1
            
            disp('steady-state inertia')
            
            % Body kinematics in strokeplane reference frame:
            
            v_st = zeros(3,n);
            a_st = zeros(3,n); 
            w_st = zeros(3,n);
            w_dot_st = zeros(3,n);
            
            [ Inertia_f ] = Inertia_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st);
            
            F_acc_b = Inertia_f.F_acc_b;
            F_vel_b = Inertia_f.F_vel_b;
            M_acc_b = Inertia_f.M_acc_b;
            M_vel_b = Inertia_f.M_vel_b;
            F_acc_st = Inertia_f.F_acc_st;
            F_vel_st = Inertia_f.F_vel_st;
            M_acc_st = Inertia_f.M_acc_st;
            M_vel_st = Inertia_f.M_vel_st;
            F_I_b = Inertia_f.F_I_b;
            M_I_b = Inertia_f.M_I_b;
            F_I_st = Inertia_f.F_I_st;
            M_I_st = Inertia_f.M_I_st;
            Lin_momentum_b = Inertia_f.Lin_momentum_b;
            Ang_momentum_b = Inertia_f.Ang_momentum_b;
            Lin_momentum_st = Inertia_f.Lin_momentum_st;
            Ang_momentum_st = Inertia_f.Ang_momentum_st;
            
            figure()
            subplot(2,1,1); plot(time,F_I_b(1,:),time,F_I_b(2,:),time,F_I_b(3,:))
            title('Inertial force in body reference frame')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_x','F_y','F_z')
            subplot(2,1,2); plot(time,F_I_st(1,:),time,F_I_st(2,:),time,F_I_st(3,:))
            title('Inertial force in strokeplane reference frame')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_x','F_y','F_z')
            
            figure()
            subplot(2,1,1); plot(time,M_I_b(1,:),time,M_I_b(2,:),time,M_I_b(3,:))
            title('Inertial moment in body reference frame')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_x','M_y','M_z')
            subplot(2,1,2); plot(time,M_I_st(1,:),time,M_I_st(2,:),time,M_I_st(3,:))
            title('Inertial moment in strokeplane reference frame')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_x','M_y','M_z')
            
            figure()
            subplot(3,1,1); plot(time,F_acc_b(1,:),time,F_vel_b(1,:),time,F_I_b(1,:))
            title('Inertial moment in x-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,2); plot(time,F_acc_b(2,:),time,F_vel_b(2,:),time,F_I_b(2,:))
            title('Inertial moment in y-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,3); plot(time,F_acc_b(3,:),time,F_vel_b(3,:),time,F_I_b(3,:))
            title('Inertial moment in z-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            
            
        case 2
            
            disp('w_x_st = 50 rad/s')
            
            % Body kinematics in strokeplane reference frame:
            
            v_st = zeros(3,n);
            a_st = zeros(3,n); 
            w_st = zeros(3,n);
            
            w_st(1,:) = 50;
            
            w_dot_st = zeros(3,n);
            
            [ Inertia_f ] = Inertia_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st);
            
            F_acc_b = Inertia_f.F_acc_b;
            F_vel_b = Inertia_f.F_vel_b;
            M_acc_b = Inertia_f.M_acc_b;
            M_vel_b = Inertia_f.M_vel_b;
            F_acc_st = Inertia_f.F_acc_st;
            F_vel_st = Inertia_f.F_vel_st;
            M_acc_st = Inertia_f.M_acc_st;
            M_vel_st = Inertia_f.M_vel_st;
            F_I_b = Inertia_f.F_I_b;
            M_I_b = Inertia_f.M_I_b;
            F_I_st = Inertia_f.F_I_st;
            M_I_st = Inertia_f.M_I_st;
            Lin_momentum_b = Inertia_f.Lin_momentum_b;
            Ang_momentum_b = Inertia_f.Ang_momentum_b;
            Lin_momentum_st = Inertia_f.Lin_momentum_st;
            Ang_momentum_st = Inertia_f.Ang_momentum_st;
            
            figure()
            subplot(2,1,1); plot(time,F_I_b(1,:),time,F_I_b(2,:),time,F_I_b(3,:))
            title('Inertial force in body reference frame')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_x','F_y','F_z')
            subplot(2,1,2); plot(time,F_I_st(1,:),time,F_I_st(2,:),time,F_I_st(3,:))
            title('Inertial force in strokeplane reference frame')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_x','F_y','F_z')
            
            figure()
            subplot(2,1,1); plot(time,M_I_b(1,:),time,M_I_b(2,:),time,M_I_b(3,:))
            title('Inertial moment in body reference frame')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_x','M_y','M_z')
            subplot(2,1,2); plot(time,M_I_st(1,:),time,M_I_st(2,:),time,M_I_st(3,:))
            title('Inertial moment in strokeplane reference frame')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_x','M_y','M_z')
            
            figure()
            subplot(3,1,1); plot(time,F_acc_b(1,:),time,F_vel_b(1,:),time,F_I_b(1,:))
            title('Inertial force in x-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,2); plot(time,F_acc_b(2,:),time,F_vel_b(2,:),time,F_I_b(2,:))
            title('Inertial force in y-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,3); plot(time,F_acc_b(3,:),time,F_vel_b(3,:),time,F_I_b(3,:))
            title('Inertial force in z-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')         
            
            figure()
            subplot(3,1,1); plot(time,M_acc_b(1,:),time,M_vel_b(1,:),time,M_I_b(1,:))
            title('Inertial moment in x-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')
            subplot(3,1,2); plot(time,M_acc_b(2,:),time,M_vel_b(2,:),time,M_I_b(2,:))
            title('Inertial moment in y-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')
            subplot(3,1,3); plot(time,M_acc_b(3,:),time,M_vel_b(3,:),time,M_I_b(3,:))
            title('Inertial moment in z-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')  
            
            
        case 3
            
            disp('w_dot_x_st = 4e4 rad/s^2 & w_x_st = 100 rad/s')
            
            % Body kinematics in strokeplane reference frame:
            
            v_st = zeros(3,n);
            
            v_st(1,:) = 200;
            
            a_st = zeros(3,n);
            
            w_st = zeros(3,n);
            
%             w_st(1,:) = 100;
            
            w_dot_st = zeros(3,n);
            
%             w_dot_st(3,:) = 2e4;
            
            
            [ Inertia_f ] = Inertia_forces( wing_kin, Body, n, v_st, a_st, w_st, w_dot_st);
            
            F_acc_b = Inertia_f.F_acc_b;
            F_vel_b = Inertia_f.F_vel_b;
            M_acc_b = Inertia_f.M_acc_b;
            M_vel_b = Inertia_f.M_vel_b;
            F_acc_st = Inertia_f.F_acc_st;
            F_vel_st = Inertia_f.F_vel_st;
            M_acc_st = Inertia_f.M_acc_st;
            M_vel_st = Inertia_f.M_vel_st;
            F_I_b = Inertia_f.F_I_b;
            M_I_b = Inertia_f.M_I_b;
            F_I_st = Inertia_f.F_I_st;
            M_I_st = Inertia_f.M_I_st;
            Lin_momentum_b = Inertia_f.Lin_momentum_b;
            Ang_momentum_b = Inertia_f.Ang_momentum_b;
            Lin_momentum_st = Inertia_f.Lin_momentum_st;
            Ang_momentum_st = Inertia_f.Ang_momentum_st;

            
            figure()
            subplot(2,1,1); plot(time,F_I_b(1,:),time,F_I_b(2,:),time,F_I_b(3,:))
            title('Inertial force in body reference frame')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_x','F_y','F_z')
            subplot(2,1,2); plot(time,F_I_st(1,:),time,F_I_st(2,:),time,F_I_st(3,:))
            title('Inertial force in strokeplane reference frame')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_x','F_y','F_z')
            
            figure()
            subplot(2,1,1); plot(time,M_I_b(1,:),time,M_I_b(2,:),time,M_I_b(3,:))
            title('Inertial moment in body reference frame')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_x','M_y','M_z')
            subplot(2,1,2); plot(time,M_I_st(1,:),time,M_I_st(2,:),time,M_I_st(3,:))
            title('Inertial moment in strokeplane reference frame')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_x','M_y','M_z')
            
            figure()
            subplot(3,1,1); plot(time,F_acc_b(1,:),time,F_vel_b(1,:),time,F_I_b(1,:))
            title('Inertial force in x-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,2); plot(time,F_acc_b(2,:),time,F_vel_b(2,:),time,F_I_b(2,:))
            title('Inertial force in y-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,3); plot(time,F_acc_b(3,:),time,F_vel_b(3,:),time,F_I_b(3,:))
            title('Inertial force in z-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')         
            
            figure()
            subplot(3,1,1); plot(time,M_acc_b(1,:),time,M_vel_b(1,:),time,M_I_b(1,:))
            title('Inertial moment in x-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')
            subplot(3,1,2); plot(time,M_acc_b(2,:),time,M_vel_b(2,:),time,M_I_b(2,:))
            title('Inertial moment in y-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')
            subplot(3,1,3); plot(time,M_acc_b(3,:),time,M_vel_b(3,:),time,M_I_b(3,:))
            title('Inertial moment in z-direction, body reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')    
            
            figure()
            subplot(3,1,1); plot(time,F_acc_st(1,:),time,F_vel_st(1,:),time,F_I_st(1,:))
            title('Inertial force in x-direction, strokeplane reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,2); plot(time,F_acc_st(2,:),time,F_vel_st(2,:),time,F_I_st(2,:))
            title('Inertial force in y-direction, strokeplane reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')
            subplot(3,1,3); plot(time,F_acc_st(3,:),time,F_vel_st(3,:),time,F_I_st(3,:))
            title('Inertial force in z-direction, strokeplane reference frame.')
            xlabel('time [s]')
            ylabel('Force in [mN]')
            legend('F_{acc}','F_{vel}','F_{tot}')         
            
            figure()
            subplot(3,1,1); plot(time,M_acc_st(1,:),time,M_vel_st(1,:),time,M_I_st(1,:),time,ones(n,1).*mean(M_vel_st(1,:)))
            title('Inertial moment in x-direction, strokeplane reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')
            subplot(3,1,2); plot(time,M_acc_st(2,:),time,M_vel_st(2,:),time,M_I_st(2,:),time,ones(n,1).*mean(M_vel_st(2,:)))
            title('Inertial moment in y-direction, strokeplane reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}')
            subplot(3,1,3); plot(time,M_acc_st(3,:),time,M_vel_st(3,:),time,M_I_st(3,:),time,ones(n,1).*mean(M_vel_st(3,:)))
            title('Inertial moment in z-direction, strokeplane reference frame.')
            xlabel('time [s]')
            ylabel('Moment in [mN*mm]')
            legend('M_{acc}','M_{vel}','M_{tot}') 
            
            
        case 4
            
            disp('conservation of momentum test')
            
            
            
        otherwise
            
            disp('no case selected')
            
    end
            


end

