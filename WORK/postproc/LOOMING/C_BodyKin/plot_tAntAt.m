% 1st order polyfit

t_AnAt = [t_Anresp t_Atresp];
t_AnAt = sortrows(t_AnAt);

nonan = isnan(t_AnAt(:,1)) + isnan(t_AnAt(:,2));
t_AnAt = t_AnAt(nonan==0,:);

t_An_sort = t_AnAt(:,1);
t_At_sort = t_AnAt(:,2);

[p_t,S_t] = polyfit(t_An_sort,t_At_sort,1);
t_At_sort_polyfit = polyval(p_t,t_An_sort);

p=0
dn=10
dm=5
clear dt_An_std_now dt_An_ste_now t_An_now t_At_now t_At_high t_At_low dt_At_ste_now
for n = [1:dm:(length(t_At_sort)-dn) length(t_At_sort)-dn]
    
    t_An_sub = t_An_sort(n:n+dn-1);
    t_At_sub = t_At_sort(n:n+dn-1);
    
    t_At_polyfit_sub = t_At_sort_polyfit(n:n+dn-1);
    dt_At_sub = t_At_sub - t_At_polyfit_sub;

    p=p+1;
    dt_At_std_now(p) = nanstd(dt_At_sub);
    dt_At_ste_now(p) = nanstd(dt_At_sub)./sqrt(length(dt_At_sub));
    
    t_An_now(p) = nanmean(t_An_sub);
    t_At_now(p) = nanmean(t_At_sub);
end

t_At_high = t_At_now + 1.96*dt_At_ste_now;
t_At_low = t_At_now - 1.96*dt_At_ste_now;

% plot(0,0,'.r','markersize',25)
% hold on
% plot(0,0,'.g','markersize',25)
% plot(0,0,'.b','markersize',25)
% legend('165deg','64deg','64deg OFF') 

% hold off
ciplot(t_At_low,t_At_high,t_An_now,[.5 .5 .5])
hold on
alpha(.25)
plot(t_An_now,t_At_now,'-k')
plot([0,max(t_AnAt(:))],[0,max(t_AnAt(:))],'--k')

% plot(t_Anresp,t_Atresp,'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',5)
plot(t_Anresp,t_Atresp,'ok','MarkerFaceColor',plotcolor,'MarkerSize',5)
% for i=1:size(t_Anresp,2)
%     if  settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 165
%         plot(t_Atresp(i),t_Anresp(i),'.r','markersize',25)
%     elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 0
%         plot(t_Atresp(i),t_Anresp(i),'.g','markersize',25)
%     elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 1
%         plot(t_Atresp(i),t_Anresp(i),'.b','markersize',25)
%     end
% end
% grid on
% axis equal

