%% plot flightpath timeline tshift headingstart

% close all
figure
subplot(3,2,1)
hold on
subplot(3,2,3)
hold on
subplot(3,2,5)
hold on
subplot(3,2,2)
hold on
subplot(3,2,4)
hold on
subplot(3,2,6)
hold on

for j=1:length(subset_seqs)
    i=subset_seqs(j);
    size(stim_angle_vel,2) - i;
    
    if isnan(color_var(i)) == 0 && color_var(i)~=0
        subplot(3,2,1)
        plot(t-t_shift(i),stim_angle_vel_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,3)
        plot(t-t_shift(i),dV_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,5)
        plot(t-t_shift(i),stim_angle_accel_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        
        subplot(3,2,2)
        plot(t-t_shift(i),F_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,4)
        plot(t-t_shift(i),Fsp_roll_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
        subplot(3,2,6)
        plot(t-t_shift(i),Fsp_pitch_plot(:,i),'-','color',cmap_plot(color_var(i),:),'linewidth',linewidth_timelines)
    end
end

subplot(3,2,1)
axis([t_start t_stop -180 180])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-180:90:180,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('heading [deg]','fontsize',10) 
%    grid on 
subplot(3,2,3)
axis([t_start t_stop -.25 .5])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',-.25:.25:1,'fontsize',8)
%     xlabel('time','fontsize',10)
ylabel('dV [m/s]','fontsize',10) 
%    grid on 
subplot(3,2,5)
axis([t_start t_stop -180 180])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',-180:90:180,'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('dir F_h_o_r [deg]','fontsize',10)
%    grid on 

subplot(3,2,2)
axis([t_start t_stop 0.5 2])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[0:1:2],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F/Mg [-]','fontsize',10)
%    grid on 
subplot(3,2,4)
axis([t_start t_stop -30 30])
%     set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'XTick',-0.5:.025:.5,'XTickLabel',[]) 
    set(gca,'YTick',[-30:30:30],'fontsize',8)
%     xlabel('time','fontsize',10)
    ylabel('F/Mg_r_o_l_l [deg]','fontsize',10)
%    grid on 
subplot(3,2,6)
axis([t_start t_stop -30 30])
    set(gca,'XTick',-0.5:.025:.5) 
    set(gca,'YTick',[-30:30:30],'fontsize',8)
    xlabel('time','fontsize',10)
    ylabel('F/Mg_p_i_t_c_h [deg]','fontsize',10)
%    grid on 

