% make borf plots F&M timelines

% clear
clc
close all

spanloc = .7;
Fscale = 10;
Fref = 1;
dref = 10;
c_plot = 20;

Nsteps = 100;
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
% load(FDB)

% mod_value = 0;
% F_color = [0 0 0];
% cord_color = [.5 .5 .5];

mod_value = max(FenhMods);
F_color = [1 0 0];
cord_color = [.5 .5 .5];

mod_values = [min(FenhMods) max(FenhMods)];
F_map = [.5 .5 .5; 1 0 0];
cord_map = [.75 .75 .75;0 0 0];


figure
mkdir(['Fenh_mov_figs'])
cd(['Fenh_mov_figs'])
        
for j = 1:length(steps)

    hold off
    for i = 1:length(mod_values)

        mod_now = mod_values(i);

        mod_diff = abs(mod_value_all - mod_now);
        n = find(mod_diff==min(mod_diff));

        mod_value_now = mod_value_all(n);
        Fenh_now = mod_value_now * Fenhance_norm;

        % F ALL
        Fx_now = Fx_all(:,:,n)/Mg_fly*F_robo2fly;
        Fx_now = nanmean(Fx_now,2);
        Fx_now = Fx_now(isnan(Fx_now)==0);

        Fz_now = Fz_all(:,:,n)/Mg_fly*F_robo2fly;
        Fz_now = nanmean(Fz_now,2);
        Fz_now = Fz_now(isnan(Fz_now)==0);

        t_now = [0:1/(length(Fx_now)-1):1];
        pp = csaps(t_now,Fx_now,.9999);
        Fx_now = fnval(pp,t_now);
        pp = csaps(t_now,Fz_now,.9999);
        Fz_now = fnval(pp,t_now);

        stroke_now = -stroke_L_all(:,n);
        stroke_now = nanmean(stroke_now,2);
        stroke_now = stroke_now(isnan(stroke_now)==0);

        pitch_now = pitch_L_all(:,n);
        pitch_now = nanmean(pitch_now,2);
        pitch_now = pitch_now(isnan(pitch_now)==0);

        dev_now = dev_L_all(:,n);
        dev_now = nanmean(dev_now,2);
        dev_now = dev_now(isnan(dev_now)==0);

        N_steps = round(steps * (length(t_now)-1) +1);


        N = N_steps(j);
        
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
        
%         Fx_tip = x_wing - Fscale*Fx_now(N);
%         Fy_tip = y_wing - Fscale*Fz_now(N);
        Fx_tip = -Fscale*Fx_now(N);
        Fy_tip = -Fscale*Fz_now(N);
        
%         plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
%         plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
%         plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
%         subplot(2,1,1)
%         hold off
%         plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_map(i,:))
        quiver(x_wing,y_wing,Fx_tip,Fy_tip,'linewidth',1,'color',F_map(i,:))
        hold on
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_map(i,:))
        plot([x_LE],[y_LE],'o','markersize',5,'markerfacecolor',cord_map(i,:),'color',cord_map(i,:))
        axis tight
        axis equal
        axis([-90 90 -30 30]) 
        axis off
        
    end
    saveas(gca, ['profNforce_Fenhance_',num2str(j),'.png'])
end
cd ..
