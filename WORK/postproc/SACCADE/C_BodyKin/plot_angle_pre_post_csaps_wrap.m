% calc angle_post
angle_pre_calc = angle_pre;
angle_post_calc = angle_post;

angle_pre_ext = angle_pre_calc;
angle_post_ext = angle_post_calc;

% wrap around
angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc -360;
angle_post_ext(end+1:end+length(angle_post_calc)) = -angle_post_calc;

angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc +360;
angle_post_ext(end+1:end+length(angle_post_calc)) = -angle_post_calc;

% average
angle_preNangle_post = [angle_pre_ext angle_post_ext];

% remove nans
angle_preNangle_post_nans = [angle_pre_ext.*angle_post_ext];
angle_preNangle_post(isnan(angle_preNangle_post_nans)==1,:) = [];

angle_preNangle_post = sortrows(angle_preNangle_post);

% [p_angle_post,S_angle_post] = polyfit(angle_preNangle_post(:,1),angle_preNangle_post(:,2),order);
pp = csaps(angle_preNangle_post(:,1),angle_preNangle_post(:,2),csaps_filt);
angle_pre_sort = angle_preNangle_post(:,1);
angle_post_sort = angle_preNangle_post(:,2);
% angle_post_sort_polyfit = polyval(p_angle_post,angle_pre_sort);
angle_post_sort_polyfit = fnval(pp,angle_pre_sort);

% variation
p=0;
clear dangle_post_std_now dangle_post_ste_now angle_post_now angle_pre_now
for n = [1:dm:(length(angle_pre_sort)-dn) length(angle_pre_sort)-dn]
    angle_post_sub = angle_post_sort(n:n+dn-1);
    angle_post_polyfit_sub = angle_post_sort_polyfit(n:n+dn-1);
    angle_pre_sub = angle_pre_sort(n:n+dn-1);
    dangle_post_sub = angle_post_sub - angle_post_polyfit_sub;

    p=p+1;
    dangle_post_std_now(p) = nanstd(dangle_post_sub);
    dangle_post_ste_now(p) = nanstd(dangle_post_sub)./sqrt(length(dangle_post_sub));
    
%     angle_post_now(p) = mean(angle_post_sub); % se not from polyfit
    angle_post_now(p) = mean(angle_post_polyfit_sub);
    
    angle_pre_now(p) = mean(angle_pre_sub);
end

% % 1st and last datapoint at edges
% angle_pre_now(1) = min(angle_pre_sort);
% angle_pre_now(end) = max(angle_pre_sort);

% 95% conf int
angle_post_high = angle_post_now + 1.96*dangle_post_ste_now;
angle_post_low = angle_post_now - 1.96*dangle_post_ste_now;

% hold off
ciplot(angle_post_low',angle_post_high',angle_pre_now',plotcolor)
hold on
alpha(.25)
% plot(angle_pre_now,angle_post_now,'-','color',plotcolor)
fnplt(pp,'color',plotcolor)

plot(angle_pre,angle_post,'ok','MarkerFaceColor',plotcolor,'markersize',5)
% plot(angle_pre_ext,angle_post_ext,'.k','MarkerFaceColor',plotcolor,'markersize',5)
axis equal

