% make borf dataset F_roll_LR

clear
clc
close all


%% NORM data
Nm2mNmm = 1e6;

%% colormap
imax = 21;

% black2blue
colormap_black2blue(1:imax,1) = 0;
colormap_black2blue(1:imax,2) = [0:.5/(imax-1):.5];
colormap_black2blue(1:imax,3) = [0:1/(imax-1):1];

%% load norm data
load('norm_data.mat')

%% Fenhance

FDB = dir('borf_db_Fenh*')
FDB = FDB.name
load(FDB)

colormap_black2blue(1:imax,1) = 0;
colormap_black2blue(1:imax,2) = [0:.5/(imax-1):.5];
colormap_black2blue(1:imax,3) = [0:1/(imax-1):1];

figure
subplot(3,3,1)
hold on
plot(Fenhance_norm*mod_value_all+1,F_mean_all/F_steady,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_freq+1,F_mean_freq/F_steady,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_stroke+1,F_mean_stroke/F_steady,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_pitch+1,F_mean_pitch/F_steady,'-k','linewidth',1)
plot(Fenhance_norm*mod_value_dev+1,F_mean_dev/F_steady,'-k','linewidth',1)
% plot(Fenhance_norm*mod_value_all(1:MODmin),F_mean_sum/F_steady,'-k','linewidth',1)

% all
for i = 2:2:length(mod_value_all)
    mod_now = mod_value_all(i);
    if mod_now>=0 && mod_now<=2
        plot(Fenhance_norm*mod_now+1,F_mean_all(i)/F_steady,'o',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% freq
for i = 1:2:length(mod_value_freq)
    mod_now = mod_value_freq(i);
    if mod_now>=0 && mod_now<=2
        plot(Fenhance_norm*mod_now+1,F_mean_freq(i)/F_steady,'v',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% stroke
for i = 1:2:length(mod_value_stroke)
    mod_now = mod_value_stroke(i);
    if mod_now>=0 && mod_now<=2
        plot(Fenhance_norm*mod_now+1,F_mean_stroke(i)/F_steady,'s',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% pitch
for i = 1:2:length(mod_value_pitch)
    mod_now = mod_value_pitch(i);
    if mod_now>=0 && mod_now<=2
        plot(Fenhance_norm*mod_now+1,F_mean_pitch(i)/F_steady,'^',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% dev
for i = 2:2:length(mod_value_dev)
    mod_now = mod_value_dev(i);
    if mod_now>=0 && mod_now<=2
        plot(Fenhance_norm*mod_now+1,F_mean_dev(i)/F_steady,'d',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

xlabel('F/Mg in (mNmm)','fontsize',10) 
ylabel('F/Mg out (mNmm)','fontsize',10)

x_min = .8;
x_max = 2;
dx = .2;
y_min = .8;
y_max = 2;
dy = .2;

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[1:dx:x_max],'fontsize',8) 
set(gca,'YTick',[1:dy:y_max],'fontsize',8) 


%% Mroll

rollDB = dir('borf_db_F_roll*')
rollDB = rollDB.name
load(rollDB)

subplot(3,3,4)
hold on
plot(Mroll_norm*mod_value_all,Nm2mNmm*(Mx_mean_all-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_freq,Nm2mNmm*(Mx_mean_freq-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_stroke,Nm2mNmm*(Mx_mean_stroke-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_pitch,Nm2mNmm*(Mx_mean_pitch-Mx_steady),'-k','linewidth',1)
plot(Mroll_norm*mod_value_dev,Nm2mNmm*(Mx_mean_dev-Mx_steady),'-k','linewidth',1)
% plot(Mroll_norm*mod_value_all(1:MODmin),Nm2mNmm*(Mx_mean_sum-Mx_steady),'-k','linewidth',1)

% all
for i = 2:2:length(mod_value_all)
    mod_now = mod_value_all(i);
    if mod_now>=0 && mod_now<=2
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_all(i)-Mx_steady),'o',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% freq
for i = 2:2:length(mod_value_freq)
    mod_now = mod_value_freq(i);
    if mod_now>=0 && mod_now<=2
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_freq(i)-Mx_steady),'v',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% stroke
for i = 1:2:length(mod_value_stroke)
    mod_now = mod_value_stroke(i);
    if mod_now>=0 && mod_now<=2
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_stroke(i)-Mx_steady),'s',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% pitch
for i = 2:2:length(mod_value_pitch)
    mod_now = mod_value_pitch(i);
    if mod_now>=0 && mod_now<=2
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_pitch(i)-Mx_steady),'^',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% dev
for i = 1:2:length(mod_value_dev)
    mod_now = mod_value_dev(i);
    if mod_now>=0 && mod_now<=2
        plot(Mroll_norm*mod_now,Nm2mNmm*(Mx_mean_dev(i)-Mx_steady),'d',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

xlabel('Mroll in (mNmm)','fontsize',10) 
ylabel('Mroll out (mNmm)','fontsize',10)

x_min = -.5e-3;
x_max = 6e-3;
y_min = -.5e-3;
y_max = 6e-3;
dy = 1e-3;

axis equal
axis([x_min x_max y_min y_max])
set(gca,'XTick',[0:dy:x_max],'fontsize',8) 
set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

%% Mpitch

pitchDB = dir('borf_db_Pitch*')
pitchDB = pitchDB.name
load(pitchDB)

subplot(3,3,7)
hold on
plot(Mpitch_norm*mod_value_all,Nm2mNmm*(My_mean_all-My_steady),'-k','linewidth',1)
plot(Mpitch_norm*mod_value_freq,Nm2mNmm*(My_mean_freq-My_steady),'-k','linewidth',1)
plot(Mpitch_norm*mod_value_stroke,Nm2mNmm*(My_mean_stroke-My_steady),'-k','linewidth',1)
plot(Mpitch_norm*mod_value_pitch,Nm2mNmm*(My_mean_pitch-My_steady),'-k','linewidth',1)
plot(Mpitch_norm*mod_value_dev,Nm2mNmm*(My_mean_dev-My_steady),'-k','linewidth',1)
% plot(Mpitch_norm*mod_value_all(1:MODmin),Nm2mNmm*(My_mean_sum-My_steady),'-k','linewidth',1)

% all
for i = 2:2:length(mod_value_all)
    mod_now = mod_value_all(i);
    if mod_now>=0 && mod_now<=2
        plot(Mpitch_norm*mod_now,Nm2mNmm*(My_mean_all(i)-My_steady),'o',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% freq
for i = 2:2:length(mod_value_freq)
    mod_now = mod_value_freq(i);
    if mod_now>=0 && mod_now<=2
        plot(Mpitch_norm*mod_now,Nm2mNmm*(My_mean_freq(i)-My_steady),'v',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% stroke
for i = 1:2:length(mod_value_stroke)
    mod_now = mod_value_stroke(i);
    if mod_now>=0 && mod_now<=2
        plot(Mpitch_norm*mod_now,Nm2mNmm*(My_mean_stroke(i)-My_steady),'s',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% pitch
for i = 1:2:length(mod_value_pitch)
    mod_now = mod_value_pitch(i);
    if mod_now>=0 && mod_now<=2
        plot(Mpitch_norm*mod_now,Nm2mNmm*(My_mean_pitch(i)-My_steady),'^',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

% dev
for i = 2:2:length(mod_value_dev)
    mod_now = mod_value_dev(i);
    if mod_now>=0 && mod_now<=2
        plot(Mpitch_norm*mod_now,Nm2mNmm*(My_mean_dev(i)-My_steady),'d',...
            'markerfacecolor',colormap_black2blue(round(10*mod_now)+1,:),...
            'markeredgecolor','k','markersize',4)
    end
end

xlabel('Mpitch in (mNmm)','fontsize',10) 
ylabel('Mpitch out (mNmm)','fontsize',10)

x_min = -.5e-3;
x_max = 4e-3;
y_min = -.5e-3;
y_max = 4e-3;
dy = 1e-3;

% axis equal
% axis([x_min x_max y_min y_max])
% set(gca,'XTick',[0:dy:x_max],'fontsize',8) 
% set(gca,'YTick',[0:dy:y_max],'fontsize',8) 

axis equal
axis([-x_max x_max -y_max y_max])
set(gca,'XTick',[-x_max:dy:x_max],'fontsize',8) 
set(gca,'YTick',[-x_max:dy:y_max],'fontsize',8) 

%% save plot

saveas(gca, 'MSfig_Fenhance_Mroll_Mpitch.fig')
saveas(gca, 'MSfig_Fenhance_Mroll_Mpitch.png')
plot2svg(['MSfig_Fenhance_Mroll_Mpitch.svg'])
    