% make borf plots F&M timelines

clear
clc
close all

spanloc = .7;
Fscale = .0002;
Fref = 1;
dref = .5e-3;

Nsteps = 8;
steps = [0:1/Nsteps:1];
steps = steps(1:end-1);

%% load norm data
load('norm_data.mat')

FenhMods = Fenhanses/Fenhance_norm
RollMods = rollaccels/rollaccel_norm
% PitchMods = pitchaccels/pitchaccel_norm
PitchMods = [-max(pitchaccels) 0 max(pitchaccels)]'/pitchaccel_norm

% FDB = 'borf_db_Fenhance_NOcali_alldata.mat'
% rollDB = 'borf_db_F_roll_LR_NOcali_means.mat'
% pitchDB = 'borf_db_PitchAccel_NOcali_means.mat'

FDB = 'borf_db_Fenhance_cali_alldata.mat'
rollDB = 'borf_db_F_roll_LR_cali_means.mat'
pitchDB = 'borf_db_PitchAccel_cali_means.mat'


%% Fenhance
load(FDB)

mod_values = FenhMods;
%%%%% !!!!!!!!!!! NEW WRONG MODS !!!!!!!!!!!!!!
% mod_values = [0 .38 .74];
% color_map = [0 0 0; 0 0 1; 0 .5 1];
mod_values = [0 .74];
F_map = [.5 .5 .5; 0 .5 1];
cord_map = [.25 .25 .25; 0 0 0];

figure
for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_all - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL
    Fx_now = Fx_all(:,:,n)/Mg_fly*F_robo2fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_all(:,:,n)/Mg_fly*F_robo2fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = stroke_L_all(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_all(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_all(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == max(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_map(i,:))
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_map(i,:))
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_map(i,:),'color',cord_map(i,:))
        
        % upstroke
        N = Nus_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_map(i,:))
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_map(i,:))
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_map(i,:),'color',cord_map(i,:))
    end
end


subplot(2,1,1)
axis equal
axis([-3e-3 3e-3 -1e-3 1e-3])

subplot(2,1,2)
plot([-2.5e-3 -2.5e-3+Fref*Fscale],[.5e-3 .5e-3],'-k','linewidth',2)
plot([+2e-3 +2e-3+dref],[.5e-3 .5e-3],'-k','linewidth',2)
axis equal
axis([-3e-3 3e-3 -1e-3 1e-3])

saveas(gca, 'MSfig_Fenhance_wingprofNforce.fig')
saveas(gca, 'MSfig_Fenhance_wingprofNforce.png')
plot2svg(['MSfig_Fenhance_wingprofNforce.svg'])


%% Mroll
load(rollDB)
figure

%% steady wing
% F_color = [.5 .5 .5];
% cord_color = [.5 .5 .5];
% mod_now = 0;
%     
%     mod_diff = abs(mod_value_allNOfreq - mod_now);
%     n = find(mod_diff==min(mod_diff));
%     
%     % F ALL NO FREQ
%     Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
%     Fx_now = nanmean(Fx_now,2);
%     Fx_now = Fx_now(isnan(Fx_now)==0);
%     
%     Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
%     Fz_now = nanmean(Fz_now,2);
%     Fz_now = Fz_now(isnan(Fz_now)==0);
%     
%     stroke_now = stroke_L_allNOfreq(:,n);
%     stroke_now = nanmean(stroke_now,2);
%     stroke_now = stroke_now(isnan(stroke_now)==0);
%     
%     pitch_now = pitch_L_allNOfreq(:,n);
%     pitch_now = nanmean(pitch_now,2);
%     pitch_now = pitch_now(isnan(pitch_now)==0);
%     
%     dev_now = dev_L_allNOfreq(:,n);
%     dev_now = nanmean(dev_now,2);
%     dev_now = dev_now(isnan(dev_now)==0);
%     
%     t_now = [0:1/(length(Fx_now)-1):1];
%     Nds_end = find(stroke_now == max(stroke_now));
%     Nus_end = length(t_now);
%     
%     Nds_steps = round(steps * (Nds_end-1) +1);
%     Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
%     
%     for j = 1:length(steps)
%         
%         % downstroke
%         N = Nds_steps(j);
%         
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
%         
%         x_LE = x_wing + dx;
%         x_TE = x_wing - dx;
%         y_LE = y_wing + dy;
%         y_TE = y_wing - dy;
%         
%         Fx_tip = x_wing - Fscale*Fx_now(N);
%         Fy_tip = y_wing - Fscale*Fz_now(N);
%         
%         subplot(2,1,1)
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
%         
%         % upstroke
%         N = Nus_steps(j);
%         
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
%         
%         x_LE = x_wing + dx;
%         x_TE = x_wing - dx;
%         y_LE = y_wing + dy;
%         y_TE = y_wing - dy;
%         
%         Fx_tip = x_wing - Fscale*Fx_now(N);
%         Fy_tip = y_wing - Fscale*Fz_now(N);
%         
%         subplot(2,1,2)
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
%     end
    
    
%% roll accel max LEFT wing
F_color = [1 0 0];
cord_color = [0 0 0];
mod_now = max(RollMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == max(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
    end
    
%% roll accel max RIGHT wing
F_color = [0 .5 1];
cord_color = [.25 .25 .25];
mod_now = max(RollMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = stroke_R_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_R_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_R_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == max(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
    end
    
subplot(2,1,1)
axis equal
axis([-3e-3 3e-3 -1e-3 1e-3])

subplot(2,1,2)
plot([-2.5e-3 -2.5e-3+Fref*Fscale],[.5e-3 .5e-3],'-k','linewidth',2)
plot([+2e-3 +2e-3+dref],[.5e-3 .5e-3],'-k','linewidth',2)
axis equal
axis([-3e-3 3e-3 -1e-3 1e-3])

saveas(gca, 'MSfig_RollAccel_wingprofNforce.fig')
saveas(gca, 'MSfig_RollAccel_wingprofNforce.png')
plot2svg(['MSfig_RollAccel_wingprofNforce.svg'])

%% Mpitch
load(pitchDB)
figure

%% steady wing
% F_color = [.5 .5 .5];
% cord_color = [.5 .5 .5];
% mod_now = 0;
%     
%     mod_diff = abs(mod_value_allNOfreq - mod_now);
%     n = find(mod_diff==min(mod_diff));
%     
%     % F ALL NO FREQ
%     Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
%     Fx_now = nanmean(Fx_now,2);
%     Fx_now = Fx_now(isnan(Fx_now)==0);
%     
%     Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
%     Fz_now = nanmean(Fz_now,2);
%     Fz_now = Fz_now(isnan(Fz_now)==0);
%     
%     stroke_now = stroke_L_allNOfreq(:,n);
%     stroke_now = nanmean(stroke_now,2);
%     stroke_now = stroke_now(isnan(stroke_now)==0);
%     
%     pitch_now = pitch_L_allNOfreq(:,n);
%     pitch_now = nanmean(pitch_now,2);
%     pitch_now = pitch_now(isnan(pitch_now)==0);
%     
%     dev_now = dev_L_allNOfreq(:,n);
%     dev_now = nanmean(dev_now,2);
%     dev_now = dev_now(isnan(dev_now)==0);
%     
%     t_now = [0:1/(length(Fx_now)-1):1];
%     Nds_end = find(stroke_now == max(stroke_now));
%     Nus_end = length(t_now);
%     
%     Nds_steps = round(steps * (Nds_end-1) +1);
%     Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
%     
%     for j = 1:length(steps)
%         
%         % downstroke
%         N = Nds_steps(j);
%         
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
%         
%         x_LE = x_wing + dx;
%         x_TE = x_wing - dx;
%         y_LE = y_wing + dy;
%         y_TE = y_wing - dy;
%         
%         Fx_tip = x_wing - Fscale*Fx_now(N);
%         Fy_tip = y_wing - Fscale*Fz_now(N);
%         
%         subplot(2,1,1)
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
%         
%         % upstroke
%         N = Nus_steps(j);
%         
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
%         
%         x_LE = x_wing + dx;
%         x_TE = x_wing - dx;
%         y_LE = y_wing + dy;
%         y_TE = y_wing - dy;
%         
%         Fx_tip = x_wing - Fscale*Fx_now(N);
%         Fy_tip = y_wing - Fscale*Fz_now(N);
%         
%         subplot(2,1,2)
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
%     end
    
    
%% pitch accel UP
F_color = [1 0 0];
cord_color = [0 0 0];
mod_now = max(PitchMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == max(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
    end
    
%% pitch accel DOWN
F_color = [0 .5 1];
cord_color = [.25 .25 .25];
mod_now = min(PitchMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly*F_robo2fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == max(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
        x_wing = spanloc*Lwing*sind(stroke_now(N));
        y_wing = spanloc*Lwing*sind(dev_now(N));
        
        dx = c_fly/2 * sind(pitch_now(N));
        dy = c_fly/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',2,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',F_color,'color',cord_color)
    end
    
subplot(2,1,1)
axis equal
axis([-3e-3 3e-3 -1e-3 1e-3])

subplot(2,1,2)
plot([-2.5e-3 -2.5e-3+Fref*Fscale],[.5e-3 .5e-3],'-k','linewidth',2)
plot([+2e-3 +2e-3+dref],[.5e-3 .5e-3],'-k','linewidth',2)
axis equal
axis([-3e-3 3e-3 -1e-3 1e-3])

saveas(gca, 'MSfig_PitchAccel_wingprofNforce.fig')
saveas(gca, 'MSfig_PitchAccel_wingprofNforce.png')
plot2svg(['MSfig_PitchAccel_wingprofNforce.svg'])

