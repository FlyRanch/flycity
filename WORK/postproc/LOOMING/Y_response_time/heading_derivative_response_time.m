
% save('heading_FhorDir_response_time_data.mat','stim_angle_vel_mirror','stim_angle_accel_mirror','t','t_shift','t_start','t_stop')
clear
clc
close all

load('heading_FhorDir_response_time_data.mat')

% y_all = stim_angle_accel_mirror;
y_all = stim_angle_vel_mirror;

t_resp_min =.01;
t_resp_max =.04;

f=188.7136;

dt = t(2)-t(1)
n_resp_min = t_resp_min/dt;
n_resp_max = t_resp_max/dt;

plot_off =0;
plot_on=1
plot_on=0

if plot_on ==1
    mkdir('heading_turn_figs')
    cd('heading_turn_figs')
end

t_resp_sub_gaus6std_all = nan(500,size(y_all,2));
y_resp_sub_gaus6std_all = nan(500,size(y_all,2));
dy_resp_sub_gaus6std_all = nan(500,size(y_all,2));

t_resp_norm_sub_gaus6std_all = nan(500,size(y_all,2));
dy_resp_norm_sub_gaus6std_all = nan(500,size(y_all,2));


for i=1:size(y_all,2)
    
        y_now = y_all(:,i);
        t_now = t-t_shift(i);
%         
%         t_sub = t_now(t_now>t_start & t_now<t_stop);
%         y_sub = y_now(t_now>t_start & t_now<t_stop);
        
        t_resp = t_now(t_now>-t_resp_min & t_now<t_resp_max);
        y_resp = y_now(t_now>-t_resp_min & t_now<t_resp_max);
        
        t_all(:,i) = t_now(i);
%         t_all_sub(1:length(t_sub),i) = t_sub;
%         y_all_sub(1:length(y_sub),i) = y_sub;
        
        dy_now = gradient(rad2deg(unwrap(deg2rad(y_now))));
        ddy_now = gradient(dy_now);
        dddy_now = gradient(ddy_now);
        
        
%         
%         dy_sub = gradient(rad2deg(unwrap(deg2rad(y_sub))));
%         ddy_sub = gradient(dy_sub);
%         dddy_sub = gradient(ddy_sub);

        dy_resp = gradient(rad2deg(unwrap(deg2rad(y_resp))));
        ddy_resp = gradient(dy_resp);
        dddy_resp = gradient(ddy_resp);

        dy_max_now = max(abs(dy_resp));
        if max(abs(dy_resp)) == max(dy_resp)
            dy_max_val = 1;
        else
            dy_max_val = -1;
        end
        
        t_dy_resp_max(i) = t_resp(abs(dy_resp) == max(abs(dy_resp)));
        t_dy_resp_max(i)
        t_now_dymax = t_now - t_dy_resp_max(i);
        t_all_dymax(1:length(t_now_dymax),i) = t_now_dymax;
        
        dy_all(1:length(dy_now),i) = dy_now;
        ddy_all(1:length(dy_now),i) = ddy_now;
        dddy_all(1:length(dy_now),i) = dddy_now;
        
%         t_sub_dymax = t_sub - t_dy_resp_max(i);
%         t_all_sub_dymax(1:length(t_sub_dymax),i) = t_sub_dymax;
        

    if plot_on ==1
        
figure(1)
            
        subplot(3,1,1)
        hold off
%         plot(t_now_dymax,y_now,'-','color','k','linewidth',.25)
        plot(t_now,y_now,'-','color','k','linewidth',.25)
        ylim([-180 180])
        grid on
        
        subplot(3,1,2)
        hold off
%         plot(t_now_dymax,dy_now,'-','color','k','linewidth',.25)
        plot(t_now,dy_now,'-','color','k','linewidth',.25)
        ylim([-max(abs(dy_now)) max(abs(dy_now))])
        grid on
        
        subplot(3,1,3)
        hold off
%         plot(t_now_dymax,ddy_now,'-','color','k','linewidth',.25)
        plot(t_now,ddy_now,'-','color','k','linewidth',.25)
        ylim([ -max(abs(ddy_now)) max(abs(ddy_now))])
        grid on
        
%         subplot(4,1,4)
%         hold off
%         plot(t_now_dymax,dddy_now,'-','color','k','linewidth',.25)
%         ylim([ -max(abs(dddy_now)) max(abs(dddy_now))])
%         grid on
    end
        
        n0 = find(t_now_dymax==0);
        if dy_max_val == 1    % POS dy_max_val
            
            % n_resp_max
            stop = 0;
            j=0;
            while stop == 0
                j = j+1;
                
                if ddy_now(n0+j) > ddy_now(n0+j-1)
                    stop = 1;
                    dn_resp_max = j;
                end
            end
            
            % n_resp_min
            stop = 0;
            j=0;
            while stop == 0
                j = j+1;
                
                if ddy_now(n0-j) < ddy_now(n0-j+1)
                    stop = 1;
                    dn_resp_min = j;
                end
            end
            
        else    % NEG dy_max_val
            
            % n_resp_max
            stop = 0;
            j=0;
            while stop == 0
                j = j+1;
                
                if ddy_now(n0+j) < ddy_now(n0+j-1)
                    stop = 1;
                    dn_resp_max = j;
                end
            end
            
            % n_resp_min
            stop = 0;
            j=0;
            while stop == 0
                j = j+1;
                
                if ddy_now(n0-j) > ddy_now(n0-j+1)
                    stop = 1;
                    dn_resp_min = j;
                end
            end
        end
        
        dn_resp_gaus2std(i) = dn_resp_max + dn_resp_min +1;
        dn_resp_gaus4std(i) = 2*dn_resp_max + 2*dn_resp_min +1;
        dn_resp_gaus6std(i) = 3*dn_resp_max + 3*dn_resp_min +1;
        
        t_resp_sub_gaus2std = t_now_dymax(n0-dn_resp_min:n0+dn_resp_max);
        y_resp_sub_gaus2std = y_now(n0-dn_resp_min:n0+dn_resp_max);
        dy_resp_sub_gaus2std = dy_now(n0-dn_resp_min:n0+dn_resp_max);
        ddy_resp_sub_gaus2std = ddy_now(n0-dn_resp_min:n0+dn_resp_max);
        dddy_resp_sub_gaus2std = dddy_now(n0-dn_resp_min:n0+dn_resp_max);
        
        t_resp_sub_gaus4std = t_now_dymax(n0-2*dn_resp_min:n0+2*dn_resp_max);
        y_resp_sub_gaus4std = y_now(n0-2*dn_resp_min:n0+2*dn_resp_max);
        dy_resp_sub_gaus4std = dy_now(n0-2*dn_resp_min:n0+2*dn_resp_max);
        ddy_resp_sub_gaus4std = ddy_now(n0-2*dn_resp_min:n0+2*dn_resp_max);
        dddy_resp_sub_gaus4std = dddy_now(n0-2*dn_resp_min:n0+2*dn_resp_max);
        
        t_resp_sub_gaus6std = t_now_dymax(n0-3*dn_resp_min:n0+3*dn_resp_max);
        y_resp_sub_gaus6std = y_now(n0-3*dn_resp_min:n0+3*dn_resp_max);
        dy_resp_sub_gaus6std = dy_now(n0-3*dn_resp_min:n0+3*dn_resp_max);
        ddy_resp_sub_gaus6std = ddy_now(n0-3*dn_resp_min:n0+3*dn_resp_max);
        dddy_resp_sub_gaus6std = dddy_now(n0-3*dn_resp_min:n0+3*dn_resp_max);
        
        t_resp_sub_gaus6std_all(1:length(t_resp_sub_gaus6std),i) = t_resp_sub_gaus6std;
        y_resp_sub_gaus6std_all(1:length(y_resp_sub_gaus6std),i) = y_resp_sub_gaus6std;
        dy_resp_sub_gaus6std_all(1:length(dy_resp_sub_gaus6std),i) = dy_resp_sub_gaus6std;
        
        t_resp_norm_sub_gaus6std_all(1:length(t_resp_sub_gaus6std),i) = 2*t_resp_sub_gaus6std/(t_resp_sub_gaus6std(end)-t_resp_sub_gaus6std(1));
        dy_resp_norm_sub_gaus6std_all(1:length(dy_resp_sub_gaus6std),i) = dy_resp_sub_gaus6std/dy_resp_sub_gaus6std(t_resp_sub_gaus6std==0);
            
    if plot_on ==1

figure(1)
        subplot(3,1,1)
        hold on
%         plot(t_resp_sub_gaus6std,y_resp_sub_gaus6std,'o')
        plot(t_resp_sub_gaus6std+t_dy_resp_max(i),y_resp_sub_gaus6std,'o')
        axis([-0.04 0.08 -180 180])
        
        subplot(3,1,2)
        hold on
%         plot(t_resp_sub_gaus6std,dy_resp_sub_gaus6std,'o')
        plot(t_resp_sub_gaus6std+t_dy_resp_max(i),dy_resp_sub_gaus6std,'o')
        axis([-0.04 0.08 -max(abs(dy_resp_sub_gaus6std)) max(abs(dy_resp_sub_gaus6std))])
        
        subplot(3,1,3)
        hold on
%         plot(t_resp_sub_gaus6std,ddy_resp_sub_gaus6std,'o')
        plot(t_resp_sub_gaus6std+t_dy_resp_max(i),ddy_resp_sub_gaus6std,'o')
        axis([-0.04 0.08 -max(abs(ddy_resp_sub_gaus6std)) max(abs(ddy_resp_sub_gaus6std))])
% 
%         subplot(4,1,4)
%         hold on
%         plot(t_resp_sub_gaus6std,dddy_resp_sub_gaus6std,'o')
%         axis([-0.02 0.02 -max(abs(dddy_resp_sub_gaus6std)) max(abs(dddy_resp_sub_gaus6std))])

saveas(gca, ['heading_turn_',num2str(i),'.png'])


%         pause
    end
end

figure(2)
subplot(1,2,1)
plot(t_resp_norm_sub_gaus6std_all,y_resp_sub_gaus6std_all,'-','color',[.5 .5 .5])
axis([-1 1 -180 180])

subplot(1,2,2)
plot(t_resp_norm_sub_gaus6std_all,dy_resp_norm_sub_gaus6std_all,'-','color',[.5 .5 .5])
axis([-1 1 -.25 1.25])

% fit gaussian to norm data
t_resp_norm_sub_gaus6std_all_list = t_resp_norm_sub_gaus6std_all(:);
t_resp_sub_gaus6std_all_list = t_resp_sub_gaus6std_all(:);
dy_resp_norm_sub_gaus6std_all_list = dy_resp_norm_sub_gaus6std_all(:);
[fitresult, gof, t_norm_gaussian_std_mean] = create_gaussian_std_fit_normdata(t_resp_norm_sub_gaus6std_all_list,dy_resp_norm_sub_gaus6std_all_list,plot_off)

hold on
plot( fitresult,'-k')

saveas(gca, ['heading_turn_all_norm.png'])

dt_resp_gaus2std = dn_resp_gaus2std * dt;
df_resp_gaus2std = dt_resp_gaus2std*f;
dt_resp_gaus4std = dn_resp_gaus4std * dt;
df_resp_gaus4std = dt_resp_gaus4std*f;
dt_resp_gaus6std = dn_resp_gaus6std * dt;
df_resp_gaus6std = dt_resp_gaus6std*f;

dt_resp_gaus2std_meanstd = mean(dt_resp_gaus2std);
df_resp_gaus2std_meanstd = mean(df_resp_gaus2std);
dt_resp_gaus4std_meanstd = mean(dt_resp_gaus4std);
df_resp_gaus4std_meanstd = mean(df_resp_gaus4std);
dt_resp_gaus6std_meanstd = mean(dt_resp_gaus6std);
df_resp_gaus6std_meanstd = mean(df_resp_gaus6std);

dt_resp_gaus2std_meanstd(2) = std(dt_resp_gaus2std)
df_resp_gaus2std_meanstd(2) = std(df_resp_gaus2std)
dt_resp_gaus4std_meanstd(2) = std(dt_resp_gaus4std)
df_resp_gaus4std_meanstd(2) = std(df_resp_gaus4std)
dt_resp_gaus6std_meanstd(2) = std(dt_resp_gaus6std)
df_resp_gaus6std_meanstd(2) = std(df_resp_gaus6std)


dt_resp_2std_meangaussian = 2*t_norm_gaussian_std_mean * dt_resp_gaus2std_meanstd;
dt_resp_4std_meangaussian = 2*t_norm_gaussian_std_mean * dt_resp_gaus4std_meanstd;
dt_resp_6std_meangaussian = 2*t_norm_gaussian_std_mean * dt_resp_gaus6std_meanstd;

df_resp_2std_meangaussian = 2*t_norm_gaussian_std_mean * df_resp_gaus2std_meanstd;
df_resp_4std_meangaussian = 2*t_norm_gaussian_std_mean * df_resp_gaus4std_meanstd;
df_resp_6std_meangaussian = 2*t_norm_gaussian_std_mean * df_resp_gaus6std_meanstd;



if plot_on ==1
cd ..
end

save('heading_response_time_results.mat')