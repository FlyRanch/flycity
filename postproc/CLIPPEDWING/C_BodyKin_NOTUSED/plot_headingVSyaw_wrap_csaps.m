% calc stim_angle_velVSyaw
angle_pre_calc = angle_pre;
stim_angle_velVSyaw_calc = stim_angle_velVSyaw;

for i=1:size(stim_angle_velVSyaw,1)
    if angle_pre(i) < -135 && stim_angle_velVSyaw(i) < -135
        stim_angle_velVSyaw_calc(i) = stim_angle_velVSyaw_calc(i) + 360;
    elseif angle_pre(i) > 135 && stim_angle_velVSyaw(i) > 135
        stim_angle_velVSyaw_calc(i) = stim_angle_velVSyaw_calc(i) - 360;
    end
end

angle_pre_ext = angle_pre_calc;
stim_angle_velVSyaw_ext = stim_angle_velVSyaw_calc;

% wrap around
angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc -360;
stim_angle_velVSyaw_ext(end+1:end+length(stim_angle_velVSyaw_calc)) = stim_angle_velVSyaw_calc;

angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc +360;
stim_angle_velVSyaw_ext(end+1:end+length(stim_angle_velVSyaw_calc)) = stim_angle_velVSyaw_calc;

% 1st order polyfit
headingNstim_angle_velVSyaw = [angle_pre_ext stim_angle_velVSyaw_ext];

% remove nans
headingNstim_angle_velVSyaw_nans = [angle_pre_ext.*stim_angle_velVSyaw_ext];
headingNstim_angle_velVSyaw(isnan(headingNstim_angle_velVSyaw_nans)==1,:) = [];

headingNstim_angle_velVSyaw = sortrows(headingNstim_angle_velVSyaw);

% [p_stim_angle_velVSyaw,S_stim_angle_velVSyaw] = polyfit(headingNstim_angle_velVSyaw(:,1),headingNstim_angle_velVSyaw(:,2),order);
pp = csaps(headingNstim_angle_velVSyaw(:,1),headingNstim_angle_velVSyaw(:,2),csaps_filt);
heading_sort = headingNstim_angle_velVSyaw(:,1);
stim_angle_velVSyaw_sort = headingNstim_angle_velVSyaw(:,2);
% stim_angle_velVSyaw_sort_polyfit = polyval(p_stim_angle_velVSyaw,heading_sort);
stim_angle_velVSyaw_sort_polyfit = fnval(pp,heading_sort);

p=0
clear dstim_angle_velVSyaw_std_now dstim_angle_velVSyaw_ste_now stim_angle_velVSyaw_now heading_now
for n = [1:dm:(length(heading_sort)-dn) length(heading_sort)-dn]
    stim_angle_velVSyaw_sub = stim_angle_velVSyaw_sort(n:n+dn-1);
    stim_angle_velVSyaw_polyfit_sub = stim_angle_velVSyaw_sort_polyfit(n:n+dn-1);
    heading_sub = heading_sort(n:n+dn-1);
    dstim_angle_velVSyaw_sub = stim_angle_velVSyaw_sub - stim_angle_velVSyaw_polyfit_sub;

    p=p+1;
    dstim_angle_velVSyaw_std_now(p) = nanstd(dstim_angle_velVSyaw_sub);
    dstim_angle_velVSyaw_ste_now(p) = nanstd(dstim_angle_velVSyaw_sub)./sqrt(length(dstim_angle_velVSyaw_sub));
    
%     stim_angle_velVSyaw_now(p) = mean(stim_angle_velVSyaw_sub); % se not from polyfit
    stim_angle_velVSyaw_now(p) = mean(stim_angle_velVSyaw_polyfit_sub);
    
    heading_now(p) = mean(heading_sub);
end

% % 1st and last datapoint at edges
% heading_now(1) = min(heading_sort);
% heading_now(end) = max(heading_sort);


stim_angle_velVSyaw_high = stim_angle_velVSyaw_now + 1.96*dstim_angle_velVSyaw_ste_now;
stim_angle_velVSyaw_low = stim_angle_velVSyaw_now - 1.96*dstim_angle_velVSyaw_ste_now;

% for i=1:size(stim_angle_velVSyaw,1)
%     if angle_pre(i) > 90 && stim_angle_velVSyaw(i) <-180
%         stim_angle_velVSyaw(i) = stim_angle_velVSyaw(i) + 360;
%     end
% end
% stim_angle_velVSyaw(stim_angle_velVSyaw>180) = stim_angle_velVSyaw(stim_angle_velVSyaw>180) - 360;
% stim_angle_velVSyaw(stim_angle_velVSyaw<-180) = stim_angle_velVSyaw(stim_angle_velVSyaw<-180) + 360;

% plot(0,0,'.r','markersize',25)
% hold on
% plot(0,0,'.g','markersize',25)
% plot(0,0,'.b','markersize',25)
% legend('165deg','64deg','64deg OFF') 

% hold off
ciplot(stim_angle_velVSyaw_low',stim_angle_velVSyaw_high',heading_now',plotcolor)
hold on
alpha(.25)
% plot(heading_now,stim_angle_velVSyaw_now,'-','color',plotcolor)
fnplt(pp,'color',plotcolor)
% plot([-180,180],[180,-180],'--k')

plot(angle_pre,stim_angle_velVSyaw,'ok','MarkerFaceColor',plotcolor,'markersize',5)
% for i=1:size(angle_pre,1)
%     if  settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 165
%         plot(angle_pre(i),stim_angle_velVSyaw(i),'.r','markersize',25)
%     elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 0
%         plot(angle_pre(i),stim_angle_velVSyaw(i),'.g','markersize',25)
%     elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 1
%         plot(angle_pre(i),stim_angle_velVSyaw(i),'.b','markersize',25)
%     end
% end
% grid on
% axis equal

