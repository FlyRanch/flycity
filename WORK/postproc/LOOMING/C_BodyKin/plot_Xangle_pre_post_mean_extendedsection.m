% calc angle_post
angle_pre_calc = angle_pre;
angle_post_calc = angle_post;

angle_pre_ext = angle_pre_calc;
angle_post_ext = angle_post_calc;

% wrap in x direction
angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc + 360;
angle_post_ext(end+1:end+length(angle_post_calc)) = angle_post_calc;

angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc - 360;
angle_post_ext(end+1:end+length(angle_post_calc)) = angle_post_calc;

% average
angle_preNangle_post = [angle_pre_ext angle_post_ext];

% remove nans
angle_preNangle_post_nans = [angle_pre_ext.*angle_post_ext];
angle_preNangle_post(isnan(angle_preNangle_post_nans)==1,:) = [];

angle_preNangle_post = sortrows(angle_preNangle_post);

angle_pre_sort = angle_preNangle_post(:,1);
angle_post_sort = angle_preNangle_post(:,2);

% circ mean & 95% conf int
p=0;
clear angle_pre_mean angle_post_mean angle_post_high angle_post_low 
for n = [1:dm:(length(angle_pre_sort)-dn) length(angle_pre_sort)-dn]
    p=p+1;
    
    angle_post_sub = deg2rad(angle_post_sort(n:n+dn-1));
    
    angle_post_mean(p,1) = nanmean(angle_post_sub);
    
    angle_post_std = nanstd(angle_post_sub);
    angle_post_high(p,1) = angle_post_mean(p,1) + 1.96*angle_post_std/sqrt(dn);
    angle_post_low(p,1) = angle_post_mean(p,1) - 1.96*angle_post_std/sqrt(dn);
    
    angle_pre_sub = angle_pre_sort(n:n+dn-1);
    angle_pre_mean(p,1) = mean(angle_pre_sub);
end

angle_post_mean = rad2deg(angle_post_mean);
angle_post_high = rad2deg(angle_post_high);
angle_post_low = rad2deg(angle_post_low);

angle_pre_conf = angle_pre_mean;
% angle_pre_conf(isnan(angle_post_high)==1)=[];
angle_post_low(isnan(angle_post_high)==1)=-180;
angle_post_high(isnan(angle_post_high)==1)=180;

% hold off
% plot(angle_pre,angle_post,'ok','MarkerFaceColor',plotcolor,'markersize',3)
plot(angle_pre_ext(angle_pre_ext>-225 & angle_pre_ext<45),angle_post_ext(angle_pre_ext>-225 & angle_pre_ext<45),'ok','MarkerFaceColor',plotcolor,'markersize',5)
hold on
ciplot(angle_post_low,angle_post_high,angle_pre_conf,plotcolor)
alpha(.25)
plot(angle_pre_mean,angle_post_mean,'-','color',plotcolor,'linewidth',2)

