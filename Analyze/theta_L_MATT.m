%% plot each seqs

% for i = seqs_plot
%     figure(i)
%     axis([0 1000 -2.5 2.5]);
%     hold on
%     plot(frames, pathDB.theta_L(:,i))
%     plot(start(:,i),zeros(length(start(:,i))),'*')
%     plot(stop(:,i),zeros(length(stop(:,i))),'.')
% %     plot(x_min_time(i,:),y_min(i,:),'*')
% end

%% Trim to wingbeats

% for i = seqs_plot % individuals
%     for j = 1:length(start(:,i)) % wingbeats
%     
%         A = start(j,i);
%         B = stop(j,i);
%               
%         if B-A > min_wingbeat % get rid of missing data
%             if A+B > 0
%                 if exist('maxframe') == 1;
%                     if B < maxframe && B > minframe % trim to before stimulus
%                        if figure_share == 1; 
%                             figure(2)
%                             hold on
%                             x = linspace(0,1,length(pathDB.theta_L(start(j,i):stop(j,i))));
%                             y = pathDB.theta_L(start(j,i):stop(j,i),i);
%                             plot(x,y,'color', color(i,1:3))
%                     
%                         elseif figure_share == 0;
%                             figure(4)
%                             hold on
%                             x = (linspace(0,1,length(pathDB.theta_L(start(j,i):stop(j,i)))));
%                             y = rad2deg(pathDB.theta_L(start(j,i):stop(j,i),i));
% 
% %                             plot(x,y,'color', color(i,1:3))
%                             plot(x,y,'color', 'blue')
%                             ylabel('stroke dev - left wing [deg]')
%                        end 
%                     end   
%                 end                         
%             end  
%         end   
%     end   
% end

%% Down and Up stroke seperated

figure(T)
hold on
axis(AX_theta);
for i = seqs_plot % individuals
    for j = 1:length(stop(:,i))-1 % wingbeats
  
        A = start(j,i);
        B = wing_min(j,i);
        C = stop(j,i);
        D = wing_min(j+1,i);

        if A > 0 && B > 0 && C > 0     
            x = linspace(0,0.5,length(pathDB.theta_L(A:B)));
            x2 = linspace(0.5,1,length(pathDB.theta_L(B:C)));
            y = rad2deg(pathDB.theta_L(A:B,i));
            y2 = rad2deg(pathDB.theta_L(B:C,i));
                  if C < maxframe && C > minframe && C-A > min_wingbeat
                         plot(x,y,'color', 'red')
                         plot(x2,y2,'color','red')
                  end

        end
        
        if C > 0 && D > 0
            
            E = linspace(0,1.5,length(pathDB.theta_L(C:D)));
            F = rad2deg(pathDB.theta_L(C:D,i));
            
            x3 = E(1:(2/3*length(E)));
            y3 = F(1:length(x3))';
            

                    if D < maxframe && D > minframe
                        plot(x3,y3,'color', 'red')
                    end
        end                
    end
end

ylabel('stroke dev - left wing [deg]')
hold off

clear A B C D E F x y x2 y2 x3 y3