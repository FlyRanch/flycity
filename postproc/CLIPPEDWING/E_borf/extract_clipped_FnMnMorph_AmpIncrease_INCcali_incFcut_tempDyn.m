
%% wing morph data            
cut_ratio(m,1) = cut_ratio_now;
cut_perc(m,1) = cut_perc_now;
cut_type(m,1) = cut_type_now;

n_total = length(t);
n_wb = n_total/Nwb;
n_start = round(wb_start*n_wb+1);
n_stop = round(wb_stop*n_wb);

% borf geometry
n_geom = find(BorfMorphCutData.cut_type == cut_type_now & BorfMorphCutData.cut_ratio == cut_ratio_now);

cut_perc_geom(m,1) = BorfMorphCutData.cut_perc(n_geom);
cut_ratio_geom(m,1) = BorfMorphCutData.cut_ratio(n_geom);
cut_type_geom(m,1) = BorfMorphCutData.cut_type(n_geom);

l_ratio(m,1) = BorfMorphCutData.WingLength_ratio(n_geom);
CoA_ratio(m,1) = BorfMorphCutData.CoA_ratio(n_geom);
A_ratio(m,1) = BorfMorphCutData.WingArea_ratio(n_geom);
S1_ratio(m,1) =BorfMorphCutData.FirstMoment_ratio(n_geom);
S2_ratio(m,1) =BorfMorphCutData.SecondMoment_ratio(n_geom);
S3_ratio(m,1) =BorfMorphCutData.ThirdMoment_ratio(n_geom);

CoA_normL(m,1) = BorfMorphCutData.CoA_norm(n_geom);
A_normL(m,1) = BorfMorphCutData.WingArea_norm(n_geom);
S1_normL(m,1) =BorfMorphCutData.FirstMoment_norm(n_geom);
S2_normL(m,1) =BorfMorphCutData.SecondMoment_norm(n_geom);
S3_normL(m,1) =BorfMorphCutData.ThirdMoment_norm(n_geom);

l_normA(m,1) = BorfMorphCutData.WingLength_normA(n_geom);
CoA_normA(m,1) = BorfMorphCutData.CoA_normA(n_geom);
S1_normA(m,1) =BorfMorphCutData.FirstMoment_normA(n_geom);
S2_normA(m,1) =BorfMorphCutData.SecondMoment_normA(n_geom);
S3_normA(m,1) =BorfMorphCutData.ThirdMoment_normA(n_geom);

        %% WB data
        if Aincrease == 1;
            stroke_cut = (radtodeg(kine(:,5)));
            dev_cut = (radtodeg(kine(:,2)));
            rot_cut = (radtodeg(kine(:,1)));
            
            stroke_intact = (radtodeg(kine(:,6)));
            dev_intact = (radtodeg(kine(:,3)));
            rot_intact = (radtodeg(kine(:,4)));
        else
            stroke_cut = (radtodeg(kine(:,6)));
            dev_cut = (radtodeg(kine(:,3)));
            rot_cut = (radtodeg(kine(:,4)));
            
            stroke_intact = (radtodeg(kine(:,5)));
            dev_intact = (radtodeg(kine(:,2)));
            rot_intact = (radtodeg(kine(:,1)));
        end
        
        AmpStroke_intact(m,1) = max(stroke_intact) - min(stroke_intact);
        AmpStroke_cut(m,1) = max(stroke_cut) - min(stroke_cut);
        Amp_ratio(m,1) = AmpStroke_cut(m,1)/AmpStroke_intact(m,1);
        
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
        time_robo = t;
        time_fly = time_robo/f_robo2fly;
        
        %% read temporal dynamics INC cali
        fx_robo =   ft(:,1) + fx_cali_L + fx_cali_R;  % thrust, fwd pos
        fy_robo =   ft(:,3) + fy_cali_L + fy_cali_R;  % sideways
        fz_robo =  -ft(:,2) + fz_cali_L + fz_cali_R; % lift, down pos

        mx_robo =   ft(:,4) + mx_cali_L + mx_cali_R;  % roll
        my_robo =   ft(:,6) + my_cali_L + my_cali_R;  % pitch, up pos
        mz_robo =  -ft(:,5) + mz_cali_L + mz_cali_R;  % yaw

        fx_fly = F_robo2fly*fx_robo;
        fy_fly = F_robo2fly*fy_robo;
        fz_fly = F_robo2fly*fz_robo;

        mx_fly = M_robo2fly*mx_robo;
        my_fly = M_robo2fly*my_robo;
        mz_fly = M_robo2fly*mz_robo;

        fx_norm = fx_fly / Mg_fly;
        fy_norm = fy_fly / Mg_fly;
        fz_norm = fz_fly / Mg_fly;

        mx_norm = mx_fly / Mg_fly / Lwing;
        my_norm = my_fly / Mg_fly / Lwing;
        mz_norm = mz_fly / Mg_fly / Lwing;
        
        %% store data
        Time_robo(:,m) = time_robo;
        Time_fly(:,m) = time_fly;
        
        Stroke_IntactWing(:,m) = stroke_intact;
        Deviation_IntactWing(:,m) = dev_intact;
        Rotation_IntactWing(:,m) = rot_intact;

        Stroke_CutWing(:,m) = stroke_cut;
        Deviation_CutWing(:,m) = dev_cut;
        Rotation_CutWing(:,m) = rot_cut;

        Fx_robo(:,m) = fx_robo;
        Fy_robo(:,m) = fy_robo;
        Fz_robo(:,m) = fz_robo;

        Mx_robo(:,m) = mx_robo;
        My_robo(:,m) = my_robo;
        Mz_robo(:,m) = mz_robo;
        
        Fx_fly(:,m) = fx_fly;
        Fy_fly(:,m) = fy_fly;
        Fz_fly(:,m) = fz_fly;

        Mx_fly(:,m) = mx_fly;
        My_fly(:,m) = my_fly;
        Mz_fly(:,m) = mz_fly;
        
        Fx_norm(:,m) = fx_norm;
        Fy_norm(:,m) = fy_norm;
        Fz_norm(:,m) = fz_norm;

        Mx_norm(:,m) = mx_norm;
        My_norm(:,m) = my_norm;
        Mz_norm(:,m) = mz_norm;
        
