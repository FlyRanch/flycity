% calc turn
angle_pre_calc = angle_pre;
turn_calc = turn;

% for i=1:size(turn,1)
%     if angle_pre(i) < -135 && turn(i) < -135
%         turn_calc(i) = turn_calc(i) + 360;
%     elseif angle_pre(i) > 135 && turn(i) > 135
%         turn_calc(i) = turn_calc(i) - 360;
%     end
% end

for i=1:size(turn,1)
    if angle_pre(i) < -90 && turn(i) < -90
        turn_calc(i) = turn_calc(i) + 360;
    elseif angle_pre(i) > 90 && turn(i) > 90
        turn_calc(i) = turn_calc(i) - 360;
    end
end

angle_pre_ext = angle_pre_calc;
turn_ext = turn_calc;

% % wrap around
% angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc -360;
% turn_ext(end+1:end+length(turn_calc)) = turn_calc +360;
% 
% angle_pre_ext(end+1:end+length(angle_pre_calc)) = angle_pre_calc +360;
% turn_ext(end+1:end+length(turn_calc)) = turn_calc -360;

% N order polyfit
headingNturn = [angle_pre_ext turn_ext];

% remove nans
headingNturn_nans = [angle_pre_ext.*turn_ext];
headingNturn(isnan(headingNturn_nans)==1,:) = [];

headingNturn = sortrows(headingNturn);

% [p_turn,S_turn] = polyfit(headingNturn(:,1),headingNturn(:,2),order);
pp = csaps(headingNturn(:,1),headingNturn(:,2),csaps_filt);
heading_sort = headingNturn(:,1);
turn_sort = headingNturn(:,2);
% turn_sort_polyfit = polyval(p_turn,heading_sort);
turn_sort_polyfit = fnval(pp,heading_sort);

p=0;
clear dturn_std_now dturn_ste_now turn_now heading_now
for n = [1:dm:(length(heading_sort)-dn) length(heading_sort)-dn]
    turn_sub = turn_sort(n:n+dn-1);
    turn_polyfit_sub = turn_sort_polyfit(n:n+dn-1);
    heading_sub = heading_sort(n:n+dn-1);
    dturn_sub = turn_sub - turn_polyfit_sub;

    p=p+1;
    dturn_std_now(p) = nanstd(dturn_sub);
    dturn_ste_now(p) = nanstd(dturn_sub)./sqrt(length(dturn_sub));
    
%     turn_now(p) = mean(turn_sub); % se not from polyfit
    turn_now(p) = mean(turn_polyfit_sub);
    
    heading_now(p) = mean(heading_sub);
end

% 1st and last datapoint at edges
heading_now(1) = min(heading_sort);
heading_now(end) = max(heading_sort);

turn_high = turn_now + 1.96*dturn_ste_now;
turn_low = turn_now - 1.96*dturn_ste_now;

% for i=1:size(turn,1)
%     if angle_pre(i) > 90 && turn(i) <-180
%         turn(i) = turn(i) + 360;
%     end
% end
turn(turn>180) = turn(turn>180) - 360;
turn(turn<-180) = turn(turn<-180) + 360;

% plot(0,0,'.r','markersize',25)
% hold on
% plot(0,0,'.g','markersize',25)
% plot(0,0,'.b','markersize',25)
% legend('165deg','64deg','64deg OFF') 

% hold off
ciplot(turn_low',turn_high',heading_now',plotcolor)
hold on
alpha(.25)
% plot(heading_now,turn_now,'-','color',plotcolor)
fnplt(pp,'color',plotcolor)
plot([-180,180],[180,-180],'--k')

plot(angle_pre,turn,'ok','MarkerFaceColor',plotcolor,'markersize',5)
% plot(angle_pre_ext,turn_ext,'.k','MarkerFaceColor',plotcolor,'markersize',5)
axis equal

