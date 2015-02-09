for i=1:size(phi_R,2)
    
%     name = [num2str(settings.seq_kin(i,1)),' seq',num2str(settings.seq_kin(i,2)),' part',num2str(settings.seq_kin(i,3))];
name = ['seq',num2str(i)];

    str_now = phi_L(:,i);
    str_now(isnan(str_now)==1)=[];
    start_now = down_time_L(:,1,i);
    start_now(isnan(start_now)==1)=[];
    stop_now = up_time_L(:,1,i);
    stop_now(isnan(stop_now)==1)=[];

    hold off
    plot(str_now)
    hold on
    plot(start_now,str_now(start_now),'o')
    plot(stop_now,str_now(stop_now),'*')
    
    str_now = phi_R(:,i);
    str_now(isnan(str_now)==1)=[];
    start_now = down_time_R(:,1,i);
    start_now(isnan(start_now)==1)=[];
    stop_now = up_time_R(:,1,i);
    stop_now(isnan(stop_now)==1)=[];

    plot(str_now,'r')
    plot(start_now,str_now(start_now),'or')
    plot(stop_now,str_now(stop_now),'*r')
    

    
%     xlim([find(isnan(str_now)==0, 1 ) find(isnan(str_now)==0, 1, 'last' )])
%     title(name)
    saveas(gca,[name,'.png'])
%     pause
end

