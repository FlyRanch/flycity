% plot each seqs

% for i = seqs_plot
%     figure(i)
%     axis([0 1000 -2.5 2.5]);
%     hold on
%     plot(frames, pathDB.phi_R(:,i),'color', color(i,1:3))
%     plot(start(:,i),zeros(length(start(:,i))),'g*')
%     plot(stop(:,i),zeros(length(stop(:,i))),'y.')
%     plot(start_min(:,i),zeros(length(start_min(:,i))),'r+')
%     plot(stop_min(:,i),zeros(length(stop_min(:,i))),'b*')
%     plot(x_min_time(i,:),y_min(i,:),'*')
% end

%% Trim to wingbeats

figure(Q)
hold on
axis(AX_phi);

for i = seqs_plot % individuals
    for j = 1:length(start(:,i)) % wingbeats
    
        A = start(j,i);
        B = stop(j,i);
              
        if B-A > min_wingbeat % get rid of missing data
            if A+B > 0
                if B < maxframe && B > minframe % trim to before stimulus
                    
                    x = linspace(0,1,length(pathDB.phi_R(start(j,i):stop(j,i))));
                    y = rad2deg(pathDB.phi_R(start(j,i):stop(j,i),i));
                    
                    x2 = linspace(0,1,length(pathDB.phi_R(stop(j,i):start(j+1,i))));
                    y2 = rad2deg(pathDB.phi_R(stop(j,i):start(j+1,i),i));
                    
%                   plot(x,y,'color', color(i,1:3))
                    plot(x,y,'color', 'blue')
                    plot(x2,y2,'color','blue')
                end
            end
        end
    end
end
    
ylabel('stroke angle - right wing [deg]')
hold off
%% Down and Up stroke seperated

% figure(1)
% hold on
% 
% for i = seqs_plot % individuals
%     for j = 1:length(start(:,i)) % wingbeats
%         
%   
%         A = start(j,i);
%         B = wing_min(j,i);
%         C = stop(j,i);
% 
%         if C-A > min_wingbeat % get rid of missing data
%             if A > 0 && B > 0 && C > 0
%                     if C < maxframe && C > minframe % trim to before stimulus
%                  
%                         x = linspace(0,0.5,length(pathDB.phi_R(A:B)));
%                         x2 = linspace(0.5,1,length(pathDB.phi_R(B:C)));
%                         y = rad2deg(pathDB.phi_R(A:B,i));
%                         y2 = rad2deg(pathDB.phi_R(B:C,i));
%                         
%                         plot(x,y,'color', 'blue')
%                         plot(x2,y2,'color','blue')
%                         ylabel('stroke angle - right wing [deg]')
%                     end
%             end
%         end
%     end
% end
% 
% hold off