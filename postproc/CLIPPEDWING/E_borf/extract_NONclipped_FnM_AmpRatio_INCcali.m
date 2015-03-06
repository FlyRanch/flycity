
        %% WB data
        stroke_cut = (radtodeg(kine(:,5)));
        stroke_intact = (radtodeg(kine(:,6)));
        
        AmpStroke_intact(m,1) = max(stroke_intact) - min(stroke_intact);
        AmpStroke_cut(m,1) = max(stroke_cut) - min(stroke_cut);
        
        Amp_ratio(m,1) = AmpStroke_cut(m,1)/AmpStroke_intact(m,1)
        
        %% cali data
        stroke_cut_round = round(radtodeg(kine(:,5)));      % right wing
        stroke_intact_round = round(radtodeg(kine(:,6)));   % left wing
        
        for n = 1:length(stroke_intact_round)
            n_L = find(stroke_intact_round(n) == cali_L.stroke_interp);
            
            fx_cali_L(n,1) = cali_L.Fx_interp(n_L);
            fy_cali_L(n,1) = cali_L.Fy_interp(n_L);
            fz_cali_L(n,1) = cali_L.Fz_interp(n_L);
        
            mx_cali_L(n,1) = cali_L.Mx_interp(n_L);
            my_cali_L(n,1) = cali_L.My_interp(n_L);
            mz_cali_L(n,1) = cali_L.Mz_interp(n_L);
        
            n_R = find(stroke_cut_round(n) == cali_R.stroke_interp);
            
            fx_cali_R(n,1) = cali_R.Fx_interp(n_R);
            fy_cali_R(n,1) = cali_R.Fy_interp(n_R);
            fz_cali_R(n,1) = cali_R.Fz_interp(n_R);
        
            mx_cali_R(n,1) = cali_R.Mx_interp(n_R);
            my_cali_R(n,1) = cali_R.My_interp(n_R);
            mz_cali_R(n,1) = cali_R.Mz_interp(n_R);
        end

        %% time        
        time_robo = t(1:n_wb);
        time(:,m) = time_robo/f_robo2fly;
        time_norm(:,m) = [0:1/(n_wb-1):1]';
        
        %% read temporal dynamics INC cali
        fx_series =   ft(:,1) + fx_cali_L + fx_cali_R;  % thrust, fwd pos
        fy_series =   ft(:,3) + fy_cali_L + fy_cali_R;  % sideways
        fz_series =  -ft(:,2) + fz_cali_L + fz_cali_R; % lift, down pos

        mx_series =   ft(:,4) + mx_cali_L + mx_cali_R;  % roll
        my_series =   ft(:,6) + my_cali_L + my_cali_R;  % pitch, up pos
        mz_series =  -ft(:,5) + mz_cali_L + mz_cali_R;  % yaw
        
        for wb = wb_start:wb_stop-1
            fx_wb(:,wb-wb_start+1) = fx_series((wb-1)*n_wb+1:wb*n_wb);
            fy_wb(:,wb-wb_start+1) = fy_series((wb-1)*n_wb+1:wb*n_wb);
            fz_wb(:,wb-wb_start+1) = fz_series((wb-1)*n_wb+1:wb*n_wb);
            
            mx_wb(:,wb-wb_start+1) = mx_series((wb-1)*n_wb+1:wb*n_wb);
            my_wb(:,wb-wb_start+1) = my_series((wb-1)*n_wb+1:wb*n_wb);
            mz_wb(:,wb-wb_start+1) = mz_series((wb-1)*n_wb+1:wb*n_wb);
        end
        fx_robo = nanmean(fx_wb,2);
        fy_robo = nanmean(fy_wb,2);
        fz_robo = nanmean(fz_wb,2);
        
        mx_robo = nanmean(mx_wb,2);
        my_robo = nanmean(my_wb,2);
        mz_robo = nanmean(mz_wb,2);

        fx_fly = F_robo2fly*fx_robo;
        fy_fly = F_robo2fly*fy_robo;
        fz_fly = F_robo2fly*fz_robo;

        mx_fly = M_robo2fly*mx_robo;
        my_fly = M_robo2fly*my_robo;
        mz_fly = M_robo2fly*mz_robo;

        fx_norm(:,m) = fx_fly / Mg_fly;
        fy_norm(:,m) = fy_fly / Mg_fly;
        fz_norm(:,m) = fz_fly / Mg_fly;

        mx_norm(:,m) = mx_fly / Mg_fly / Lwing;
        my_norm(:,m) = my_fly / Mg_fly / Lwing;
        mz_norm(:,m) = mz_fly / Mg_fly / Lwing;
        
        %% mean F&M
        Fx_norm(m,1) = mean(fx_norm(:,m));
        Fy_norm(m,1) = mean(fy_norm(:,m));
        Fz_norm(m,1) = mean(fz_norm(:,m));

        Mx_norm(m,1) = mean(mx_norm(:,m));
        My_norm(m,1) = mean(my_norm(:,m));
        Mz_norm(m,1) = mean(mz_norm(:,m));
        
        My_norm_CoG(m,1) = My_norm(m,1) - d_norm_steady*Fz_norm(m,1);
