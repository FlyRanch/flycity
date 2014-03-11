% make borf plots F&M timelines

clear
clc
close all

spanloc = .7;
Fscale = 3;
Fref = 2;
dref = 10;
c_plot = 10;

Nsteps = 10;
steps = [0:1/Nsteps:1];
steps = steps(1:end-1);

%% load MOD&norm data
load('MOD_norm_data.mat')

%% Fenhance
load(FDB_all)

mod_values = FenhMods;
F_map = [.5 .5 .5; 0 1 1];
cord_map = [.75 .75 .75; 0 0 0];

figure
subplot(2,1,1)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

subplot(2,1,2)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

for i = 1:length(mod_values)
    mod_now = mod_values(i);
    
    mod_diff = abs(mod_value_all - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL
    Fx_now = Fx_all(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_all(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_L_all(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_all(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_all(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));

        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_map(i,:))
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_map(i,:))
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_map(i,:),'color',cord_map(i,:))
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_map(i,:))
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_map(i,:))
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_map(i,:),'color',cord_map(i,:))
    end
end

mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_Fenhance_wingprofNforce_angles.fig')
saveas(gca, 'MSfig_Fenhance_wingprofNforce_angles.png')
plot2svg(['MSfig_Fenhance_wingprofNforce_angles.svg'])
cd ..

%% Mroll
load(rollDB_LnRsym)
figure
subplot(2,1,1)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

subplot(2,1,2)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

%% steady wing
% F_color = [.5 .5 .5];
% cord_color = [.75 .75 .75];
% mod_now = 0;
%     
%     mod_diff = abs(mod_value_allNOfreq - mod_now);
%     n = find(mod_diff==min(mod_diff));
%     
%     % F ALL NO FREQ
%     Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly;
%     Fx_now = nanmean(Fx_now,2);
%     Fx_now = Fx_now(isnan(Fx_now)==0);
%     
%     Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly;
%     Fz_now = nanmean(Fz_now,2);
%     Fz_now = Fz_now(isnan(Fz_now)==0);
%     
%     stroke_now = -stroke_L_allNOfreq(:,n);
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
%     Nds_end = find(stroke_now == min(stroke_now));
%     Nus_end = length(t_now);
%     
%     Nds_steps = round(steps * (Nds_end-1) +1);
%     Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
%     
%     % smoove data
%     pp = csaps(t_now,Fx_now,.9999);
%     Fx_now = fnval(pp,t_now);
%     pp = csaps(t_now,Fz_now,.9999);
%     Fz_now = fnval(pp,t_now);
% 
%     for j = 1:length(steps)
%         
%         % downstroke
%         N = Nds_steps(j);
%         
% %         x_wing = spanloc*Lwing*sind(stroke_now(N));
% %         y_wing = spanloc*Lwing*sind(dev_now(N));
% %         dx = c_fly/2 * sind(pitch_now(N));
% %         dy = c_fly/2 * cosd(pitch_now(N));
%         x_wing = stroke_now(N);
%         y_wing = dev_now(N);
%         dx = c_plot/2 * sind(pitch_now(N));
%         dy = c_plot/2 * cosd(pitch_now(N));
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
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
%         
%         % upstroke
%         N = Nus_steps(j);
%         
% %         x_wing = spanloc*Lwing*sind(stroke_now(N));
% %         y_wing = spanloc*Lwing*sind(dev_now(N));
% %         dx = c_fly/2 * sind(pitch_now(N));
% %         dy = c_fly/2 * cosd(pitch_now(N));
%         x_wing = stroke_now(N);
%         y_wing = dev_now(N);
%         dx = c_plot/2 * sind(pitch_now(N));
%         dy = c_plot/2 * cosd(pitch_now(N));
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
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
%     end
    
    
%% roll accel max LEFT wing
F_color = [0 1 1];
cord_color = [0 0 0];
mod_now = max(RollMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq_left(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq_left(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);

    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
    end
    
%% roll accel max RIGHT wing
F_color = [1 0 0];
cord_color = [.75 .75 .75];
mod_now = max(RollMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq_right(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq_right(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_R_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_R_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_R_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
    end
    
mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_RollAccel_wingprofNforce_angles.fig')
saveas(gca, 'MSfig_RollAccel_wingprofNforce_angles.png')
plot2svg(['MSfig_RollAccel_wingprofNforce_angles.svg'])
cd ..

%% Myaw
load(yawDB_LnRsym)
figure
subplot(2,1,1)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

subplot(2,1,2)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

%% steady wing
% F_color = [.5 .5 .5];
% cord_color = [.75 .75 .75];
% mod_now = 0;
%     
%     mod_diff = abs(mod_value_allNOfreq - mod_now);
%     n = find(mod_diff==min(mod_diff));
%     
%     % F ALL NO FREQ
%     Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly;
%     Fx_now = nanmean(Fx_now,2);
%     Fx_now = Fx_now(isnan(Fx_now)==0);
%     
%     Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly;
%     Fz_now = nanmean(Fz_now,2);
%     Fz_now = Fz_now(isnan(Fz_now)==0);
%     
%     stroke_now = -stroke_L_allNOfreq(:,n);
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
%     Nds_end = find(stroke_now == min(stroke_now));
%     Nus_end = length(t_now);
%     
%     Nds_steps = round(steps * (Nds_end-1) +1);
%     Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
%     
%     % smoove data
%     pp = csaps(t_now,Fx_now,.9999);
%     Fx_now = fnval(pp,t_now);
%     pp = csaps(t_now,Fz_now,.9999);
%     Fz_now = fnval(pp,t_now);
% 
%     for j = 1:length(steps)
%         
%         % downstroke
%         N = Nds_steps(j);
%         
% %         x_wing = spanloc*Lwing*sind(stroke_now(N));
% %         y_wing = spanloc*Lwing*sind(dev_now(N));
% %         dx = c_fly/2 * sind(pitch_now(N));
% %         dy = c_fly/2 * cosd(pitch_now(N));
%         x_wing = stroke_now(N);
%         y_wing = dev_now(N);
%         dx = c_plot/2 * sind(pitch_now(N));
%         dy = c_plot/2 * cosd(pitch_now(N));
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
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
%         
%         % upstroke
%         N = Nus_steps(j);
%         
% %         x_wing = spanloc*Lwing*sind(stroke_now(N));
% %         y_wing = spanloc*Lwing*sind(dev_now(N));
% %         dx = c_fly/2 * sind(pitch_now(N));
% %         dy = c_fly/2 * cosd(pitch_now(N));
%         x_wing = stroke_now(N);
%         y_wing = dev_now(N);
%         dx = c_plot/2 * sind(pitch_now(N));
%         dy = c_plot/2 * cosd(pitch_now(N));
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
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
%     end
    
    
%% yaw accel max LEFT wing
F_color = [0 1 1];
cord_color = [0 0 0];
mod_now = max(YawMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq_left(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq_left(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);

    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
    end
    
%% yaw accel max RIGHT wing
F_color = [1 0 0];
cord_color = [.75 .75 .75];
mod_now = max(YawMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq_right(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq_right(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_R_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_R_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_R_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
    end
    
mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_YawAccel_wingprofNforce_angles.fig')
saveas(gca, 'MSfig_YawAccel_wingprofNforce_angles.png')
plot2svg(['MSfig_YawAccel_wingprofNforce_angles.svg'])
cd ..

%% Mpitch
load(pitchDB_all)
figure
subplot(2,1,1)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

subplot(2,1,2)
plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -30 30]) 
axis off

%% steady wing
% F_color = [.5 .5 .5];
% cord_color = [.5 .5 .5];
% mod_now = 0;
%     
%     mod_diff = abs(mod_value_allNOfreq - mod_now);
%     n = find(mod_diff==min(mod_diff));
%     
%     % F ALL NO FREQ
%     Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly;
%     Fx_now = nanmean(Fx_now,2);
%     Fx_now = Fx_now(isnan(Fx_now)==0);
%     
%     Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly;
%     Fz_now = nanmean(Fz_now,2);
%     Fz_now = Fz_now(isnan(Fz_now)==0);
%     
%     stroke_now = -stroke_L_allNOfreq(:,n);
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
%     Nds_end = find(stroke_now == min(stroke_now));
%     Nus_end = length(t_now);
%     
%     Nds_steps = round(steps * (Nds_end-1) +1);
%     Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
%     
%     % smoove data
%     pp = csaps(t_now,Fx_now,.9999);
%     Fx_now = fnval(pp,t_now);
%     pp = csaps(t_now,Fz_now,.9999);
%     Fz_now = fnval(pp,t_now);
%     
%     for j = 1:length(steps)
%         
%         % downstroke
%         N = Nds_steps(j);
%         
% %         x_wing = spanloc*Lwing*sind(stroke_now(N));
% %         y_wing = spanloc*Lwing*sind(dev_now(N));
% %         dx = c_fly/2 * sind(pitch_now(N));
% %         dy = c_fly/2 * cosd(pitch_now(N));
%         x_wing = stroke_now(N);
%         y_wing = dev_now(N);
%         dx = c_plot/2 * sind(pitch_now(N));
%         dy = c_plot/2 * cosd(pitch_now(N));
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
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
%         
%         % upstroke
%         N = Nus_steps(j);
%         
% %         x_wing = spanloc*Lwing*sind(stroke_now(N));
% %         y_wing = spanloc*Lwing*sind(dev_now(N));
% %         dx = c_fly/2 * sind(pitch_now(N));
% %         dy = c_fly/2 * cosd(pitch_now(N));
%         x_wing = stroke_now(N);
%         y_wing = dev_now(N);
%         dx = c_plot/2 * sind(pitch_now(N));
%         dy = c_plot/2 * cosd(pitch_now(N));
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
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
%         hold on
%         plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
%         plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
%     end
    
    
%% pitch accel DOWN
F_color = [1 0 0];
cord_color = [.75 .75 .75];
mod_now = max(PitchMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
    end
    
%% pitch accel UP
F_color = [0 .5 1];
cord_color = 'k';
mod_now = min(PitchMods);
    
    mod_diff = abs(mod_value_allNOfreq - mod_now);
    n = find(mod_diff==min(mod_diff));
    
    % F ALL NO FREQ
    Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly;
    Fx_now = nanmean(Fx_now,2);
    Fx_now = Fx_now(isnan(Fx_now)==0);
    
    Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly;
    Fz_now = nanmean(Fz_now,2);
    Fz_now = Fz_now(isnan(Fz_now)==0);
    
    stroke_now = -stroke_L_allNOfreq(:,n);
    stroke_now = nanmean(stroke_now,2);
    stroke_now = stroke_now(isnan(stroke_now)==0);
    
    pitch_now = pitch_L_allNOfreq(:,n);
    pitch_now = nanmean(pitch_now,2);
    pitch_now = pitch_now(isnan(pitch_now)==0);
    
    dev_now = dev_L_allNOfreq(:,n);
    dev_now = nanmean(dev_now,2);
    dev_now = dev_now(isnan(dev_now)==0);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    Nds_end = find(stroke_now == min(stroke_now));
    Nus_end = length(t_now);
    
    Nds_steps = round(steps * (Nds_end-1) +1);
    Nus_steps = round(steps * (Nus_end-Nds_end-1) +1 +Nds_end);
    
    % smoove data
    pp = csaps(t_now,Fx_now,.9999);
    Fx_now = fnval(pp,t_now);
    pp = csaps(t_now,Fz_now,.9999);
    Fz_now = fnval(pp,t_now);
    
    for j = 1:length(steps)
        
        % downstroke
        N = Nds_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,1)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
        
        % upstroke
        N = Nus_steps(j);
        
%         x_wing = spanloc*Lwing*sind(stroke_now(N));
%         y_wing = spanloc*Lwing*sind(dev_now(N));
%         dx = c_fly/2 * sind(pitch_now(N));
%         dy = c_fly/2 * cosd(pitch_now(N));
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_plot/2 * sind(pitch_now(N));
        dy = c_plot/2 * cosd(pitch_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        subplot(2,1,2)
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_color,'color',cord_color)
    end
    
mkdir('MSfigs')
cd('MSfigs')
saveas(gca, 'MSfig_PitchAccel_wingprofNforce_angles.fig')
saveas(gca, 'MSfig_PitchAccel_wingprofNforce_angles.png')
plot2svg(['MSfig_PitchAccel_wingprofNforce_angles.svg'])
cd ..
