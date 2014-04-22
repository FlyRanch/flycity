function [ cg_b_neutral ] = cg_secant_method( body_model, wing_model, wing_kin, nr_points, dt )

    % Determine the position of the center of gravity of the body for which
    % the wingbeat averaged forces and moments will be zero, using the
    % secant method.
    
    x_cg = 0;
    
    
    R_JL = body_model.Joint_left;
    
    
    % Initial conditions:
    
    state.id        = 1;
    state.dt        = dt;
    state.vb        = [0; 0; 0];
    state.wb        = [0; 0; 0];
    state.ab        = [0; 0; 0];
    state.w_dot_b   = [0; 0; 0];
    
    [ state_new, ~ ] = RK_updater( state, body_model, wing_model, wing_kin );
    
    state.id        = 1;
    state.dt        = dt;
    state.vb        = state_new.vb_0;
    state.wb        = state_new.wb_0;
    state.ab        = state_new.ab;
    state.w_dot_b   = state_new.w_dot_b;
    
    vb      = zeros(3,nr_points);
    wb      = zeros(3,nr_points);
    ab      = zeros(3,nr_points);
    w_dot_b = zeros(3,nr_points);
    vb_0    = zeros(3,nr_points);
    wb_0    = zeros(3,nr_points);
    FI_acc  = zeros(3,nr_points);
    FI_vel  = zeros(3,nr_points);
    MI_acc  = zeros(3,nr_points);
    MI_vel  = zeros(3,nr_points);
    LM      = zeros(3,nr_points);
    AM      = zeros(3,nr_points);
       
    for j = 1:nr_points
        
        vb(:,j)      = state.vb;
        wb(:,j)      = state.wb;
        ab(:,j)      = state.ab;
        w_dot_b(:,j) = state.w_dot_b;
        
        [ state_new, FM_old ] = RK_updater( state, body_model, wing_model, wing_kin );
        
        state.id        = state_new.id;
        state.dt        = state_new.dt;
        state.vb        = state_new.vb;
        state.wb        = state_new.wb;
        state.ab        = state_new.ab;
        state.w_dot_b   = state_new.w_dot_b;
                
        FI_acc(:,j)     = FM_old.FI_acc_b;
        FI_vel(:,j)     = FM_old.FI_vel_b;
        MI_acc(:,j)     = FM_old.MI_acc_b;
        MI_vel(:,j)     = FM_old.MI_vel_b;
        
        LM(:,j)         = FM_old.LM;
        AM(:,j)         = FM_old.AM;
        
        vb_0(:,j)       = state_new.vb_0;
        wb_0(:,j)       = state_new.wb_0;
        
    end
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(FI_vel(1,:))
%     subplot(3,1,2); plot(FI_vel(2,:))
%     subplot(3,1,3); plot(FI_vel(3,:))
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); plot(MI_vel(1,:))
%     subplot(3,1,2); plot(MI_vel(2,:))
%     subplot(3,1,3); plot(MI_vel(3,:))
%     hold off
% 
%     pause

    
    FI_mean = mean(FI_vel,2);
    MI_mean = mean(MI_vel,2);
    
    % Compute x_cg(2):
    
    x_cg_comp = (MI_mean(2)-R_JL(3)*FI_mean(1))/FI_mean(3);

%     x_cg_comp = mean((MI_vel(2,:)-R_JL(3)*FI_vel(1,:))./FI_vel(3,:));
    
    x_cg_n = x_cg_comp;
    
    x_cg = [x_cg x_cg_n];
    
    % Run the secant method
    
    i = 2;
    
    while abs(x_cg(i)-x_cg(i-1)) > 0.00001
        
        i
        
        x_cg(i)
        
        body_model.cg_b = [ x_cg(i); 0; 0];
        
        state.id        = 1;
        state.dt        = dt;
        state.vb        = vb_0(:,1);
        state.wb        = wb_0(:,1);
        state.ab        = ab(:,1);
        state.w_dot_b   = w_dot_b(:,1);
        
        vb      = zeros(3,nr_points);
        wb      = zeros(3,nr_points);
        ab      = zeros(3,nr_points);
        w_dot_b = zeros(3,nr_points);
        vb_0    = zeros(3,nr_points);
        wb_0    = zeros(3,nr_points);
        FI_acc  = zeros(3,nr_points);
        FI_vel  = zeros(3,nr_points);
        MI_acc  = zeros(3,nr_points);
        MI_vel  = zeros(3,nr_points);
        LM      = zeros(3,nr_points);
        AM      = zeros(3,nr_points);
   
        % RK4-simulation:

        for j = 1:nr_points
            
            vb(:,j)      = state.vb;
            wb(:,j)      = state.wb;
            ab(:,j)      = state.ab;
            w_dot_b(:,j) = state.w_dot_b;

            [ state_new, FM_old ] = RK_updater( state, body_model, wing_model, wing_kin );

            state.id        = state_new.id;
            state.dt        = state_new.dt;
            state.vb        = state_new.vb;
            state.wb        = state_new.wb;
            state.ab        = state_new.ab;
            state.w_dot_b   = state_new.w_dot_b;

            FI_acc(:,j)     = FM_old.FI_acc_b;
            FI_vel(:,j)     = FM_old.FI_vel_b;
            MI_acc(:,j)     = FM_old.MI_acc_b;
            MI_vel(:,j)     = FM_old.MI_vel_b;

            LM(:,j)         = FM_old.LM;
            AM(:,j)         = FM_old.AM;

            vb_0(:,j)       = state_new.vb_0;
            wb_0(:,j)       = state_new.wb_0;

        end
        
        FI_mean = mean(FI_vel,2);
        MI_mean = mean(MI_vel,2);
        
%         x_cg_t = mean((MI_vel(2,:)-R_JL(3)*FI_vel(1,:))./FI_vel(3,:));
        
        if abs(FI_mean) < 1e-12 & abs(MI_mean) < 1e-12
            
            x_cg_t = x_cg_comp(i-1);
            
        else

            x_cg_t = (MI_mean(2)-R_JL(3)*FI_mean(1))/FI_mean(3);

        end
        
        x_cg_comp = [x_cg_comp x_cg_t];
        
        x_cg_n = x_cg(i)-x_cg_comp(i)*((x_cg(i-1)-x_cg(i))/(x_cg_comp(i-1)-x_cg_comp(i)));
    
        x_cg = [x_cg x_cg_n];
        
        i = i+1;
        
    end
    
    figure()
    hold on
    subplot(3,1,1); plot(vb(1,:))
    subplot(3,1,2); plot(vb(2,:))
    subplot(3,1,3); plot(vb(3,:))
    hold off
    title('vb')
    
    figure()
    hold on
    subplot(3,1,1); plot(wb(1,:))
    subplot(3,1,2); plot(wb(2,:))
    subplot(3,1,3); plot(wb(3,:))
    hold off
    title('wb')
    
    mean(vb,2)
    
    mean(wb,2)
    
    figure()
    hold on
    subplot(3,1,1); plot(ab(1,:))
    subplot(3,1,2); plot(ab(2,:))
    subplot(3,1,3); plot(ab(3,:))
    hold off
    title('ab')
    
    figure()
    hold on
    subplot(3,1,1); plot(w_dot_b(1,:))
    subplot(3,1,2); plot(w_dot_b(2,:))
    subplot(3,1,3); plot(w_dot_b(3,:))
    hold off
    title('w_dot_b')
    
    mean(ab,2)
    
    mean(w_dot_b,2)
    
        
    figure()
    hold on
    subplot(3,1,1); plot(LM(1,:))
    subplot(3,1,2); plot(LM(2,:))
    subplot(3,1,3); plot(LM(3,:))
    hold off
    title('LM')
    
    figure()
    hold on
    subplot(3,1,1); plot(AM(1,:))
    subplot(3,1,2); plot(AM(2,:))
    subplot(3,1,3); plot(AM(3,:))
    hold off
    title('AM')
    
    mean(LM,2)
    
    mean(AM,2)
    
    mean(FI_vel,2)
    
    mean(MI_vel,2)
    
    figure()
    hold on
    subplot(3,1,1); plot(FI_vel(1,:))
    subplot(3,1,2); plot(FI_vel(2,:))
    subplot(3,1,3); plot(FI_vel(3,:))
    hold off
    title('FI_vel')
    
    figure()
    hold on
    subplot(3,1,1); plot(MI_vel(1,:))
    subplot(3,1,2); plot(MI_vel(2,:))
    subplot(3,1,3); plot(MI_vel(3,:))
    hold off
    title('MI_vel')
    
    figure()
    plot(x_cg)

    pause

    cg_b_neutral = body_model.cg_b;


end