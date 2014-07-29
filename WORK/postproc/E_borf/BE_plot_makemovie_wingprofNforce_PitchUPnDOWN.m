% make borf plots F&M timelines

clear
clc
close all

spanloc = .7;
Fscale = 3;
Fref = 2;
dref = 10;
c_plot = 10;

Nsteps = 200;
steps = [0:1/Nsteps:1];
steps = steps(1:end-1);

%% load MOD&norm data
load('MOD_norm_data.mat')

%% Fenhance
load(pitchDB_all)

mod_values = [min(PitchMods) max(PitchMods)]; % NO steady
F_map = [1 0 0;0 1 1];
cord_map = [0 0 0;.5 .5 .5];

figure
mkdir(['Pitch_mov_figs'])
cd(['Pitch_mov_figs'])
        
for j = 1:length(steps)

    hold off
    for i = 1:length(mod_values)

        mod_now = mod_values(i);

        mod_diff = abs(mod_value_allNOfreq - mod_now);
        n = find(mod_diff==min(mod_diff));

        mod_value_now = mod_value_allNOfreq(n);
        Fenh_now = mod_value_now * Fenhance_norm;

        % F ALL
        Fx_now = Fx_allNOfreq(:,:,n)/Mg_fly;
        Fx_now = nanmean(Fx_now,2);
        Fx_now = Fx_now(isnan(Fx_now)==0);

        Fz_now = Fz_allNOfreq(:,:,n)/Mg_fly;
        Fz_now = nanmean(Fz_now,2);
        Fz_now = Fz_now(isnan(Fz_now)==0);

        t_now = [0:1/(length(Fx_now)-1):1];
        pp = csaps(t_now,Fx_now,.9999);
        Fx_now = fnval(pp,t_now);
        pp = csaps(t_now,Fz_now,.9999);
        Fz_now = fnval(pp,t_now);

        stroke_now = -stroke_L_allNOfreq(:,n);
        stroke_now = nanmean(stroke_now,2);
        stroke_now = stroke_now(isnan(stroke_now)==0);

        pitch_now = pitch_L_allNOfreq(:,n);
        pitch_now = nanmean(pitch_now,2);
        pitch_now = pitch_now(isnan(pitch_now)==0);

        dev_now = dev_L_allNOfreq(:,n);
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
    saveas(gca, ['profNforce_PitchAccel_',num2str(j),'.png'])
end
cd ..
