%% theta - left minus right
for i = seqs_plot % individuals
    for j = 1:length(start(:,i)) % wingbeats
    
        A = start(j,i);
        B = stop(j,i);
              
        if B-A > min_wingbeat % get rid of missing data
            if A+B > 0
                if B < maxframe && B > minframe % trim to before stimulus
                      
                    right = pathDB.theta_R(start(j,i):stop(j,i),i);
                    left = pathDB.theta_L(start(j,i):stop(j,i),i);
                    
                    delta = right-left;            
                    
                    figure(6)
                    hold on
                    x = linspace(0,1,length(delta));
                    y = delta;
                    plot(x,y,'color', color(i,1:3))
                    ylabel('theta - left minus right [deg]')
                end                     
            end
        end
       
    end
end