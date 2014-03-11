% calc headingVSyaw
angle_pre_calc = angle_pre;
headingVSyaw_calc = headingVSyaw;

for i=1:size(headingVSyaw,1)
    if angle_pre(i) < -135 && headingVSyaw(i) < -135
        headingVSyaw_calc(i) = headingVSyaw_calc(i) + 360;
    elseif angle_pre(i) > 135 && headingVSyaw(i) > 135
        headingVSyaw_calc(i) = headingVSyaw_calc(i) - 360;
    end
end

angle_pre_ext = angle_pre_calc;
headingVSyaw_ext = headingVSyaw_calc;

% % wrap around
% angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc -360;
% headingVSyaw_ext(end+1:end+length(headingVSyaw_calc)) = headingVSyaw_calc +360;
% 
% angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc +360;
% headingVSyaw_ext(end+1:end+length(headingVSyaw_calc)) = headingVSyaw_calc -360;

% 1st order polyfit
headingNheadingVSyaw = [angle_pre_ext headingVSyaw_ext];

% remove nans
headingNheadingVSyaw_nans = [angle_pre_ext.*headingVSyaw_ext];
headingNheadingVSyaw(isnan(headingNheadingVSyaw_nans)==1,:) = [];

headingNheadingVSyaw = sortrows(headingNheadingVSyaw);

% [p_headingVSyaw,S_headingVSyaw] = polyfit(headingNheadingVSyaw(:,1),headingNheadingVSyaw(:,2),order);
pp = csaps(headingNheadingVSyaw(:,1),headingNheadingVSyaw(:,2),csaps_filt);
heading_sort = headingNheadingVSyaw(:,1);
headingVSyaw_sort = headingNheadingVSyaw(:,2);
% headingVSyaw_sort_polyfit = polyval(p_headingVSyaw,heading_sort);
headingVSyaw_sort_polyfit = fnval(pp,heading_sort);

p=0
clear dheadingVSyaw_std_now dheadingVSyaw_ste_now headingVSyaw_now heading_now
for n = [1:dm:(length(heading_sort)-dn) length(heading_sort)-dn]
    headingVSyaw_sub = headingVSyaw_sort(n:n+dn-1);
    headingVSyaw_polyfit_sub = headingVSyaw_sort_polyfit(n:n+dn-1);
    heading_sub = heading_sort(n:n+dn-1);
    dheadingVSyaw_sub = headingVSyaw_sub - headingVSyaw_polyfit_sub;

    p=p+1;
    dheadingVSyaw_std_now(p) = nanstd(dheadingVSyaw_sub);
    dheadingVSyaw_ste_now(p) = nanstd(dheadingVSyaw_sub)./sqrt(length(dheadingVSyaw_sub));
    
%     headingVSyaw_now(p) = mean(headingVSyaw_sub); % se not from polyfit
    headingVSyaw_now(p) = mean(headingVSyaw_polyfit_sub);
    
    heading_now(p) = mean(heading_sub);
end

% 1st and last datapoint at edges
heading_now(1) = min(heading_sort);
heading_now(end) = max(heading_sort);


headingVSyaw_high = headingVSyaw_now + 1.96*dheadingVSyaw_ste_now;
headingVSyaw_low = headingVSyaw_now - 1.96*dheadingVSyaw_ste_now;

% for i=1:size(headingVSyaw,1)
%     if angle_pre(i) > 90 && headingVSyaw(i) <-180
%         headingVSyaw(i) = headingVSyaw(i) + 360;
%     end
% end
% headingVSyaw(headingVSyaw>180) = headingVSyaw(headingVSyaw>180) - 360;
% headingVSyaw(headingVSyaw<-180) = headingVSyaw(headingVSyaw<-180) + 360;

% plot(0,0,'.r','markersize',25)
% hold on
% plot(0,0,'.g','markersize',25)
% plot(0,0,'.b','markersize',25)
% legend('165deg','64deg','64deg OFF') 

% hold off
ciplot(headingVSyaw_low',headingVSyaw_high',heading_now',plotcolor)
hold on
alpha(.25)
% plot(heading_now,headingVSyaw_now,'-','color',plotcolor)
fnplt(pp,'color',plotcolor)
% plot([-180,180],[180,-180],'--k')

plot(angle_pre,headingVSyaw,'ok','MarkerFaceColor',plotcolor,'markersize',5)
% for i=1:size(angle_pre,1)
%     if  settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 165
%         plot(angle_pre(i),headingVSyaw(i),'.r','markersize',25)
%     elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 0
%         plot(angle_pre(i),headingVSyaw(i),'.g','markersize',25)
%     elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 1
%         plot(angle_pre(i),headingVSyaw(i),'.b','markersize',25)
%     end
% end
% grid on
% axis equal

