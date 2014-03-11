%% phi - left minus right

minframe = 300;
maxframe = 1000;

fig = figure(31);
hold on
axis([0 1 0 30]);

for i = seqs_plot % individuals
    for j = 1:length(start(:,i)) % wingbeats
    
        A = start(j,i);
        B = stop(j,i);
              
        if B-A > min_wingbeat % get rid of missing data
            if A+B > 0
                if B < maxframe && B > minframe % trim to before stimulus
                if B < trims(i,3) && B > trims(i,2)
                    right = pathDB.phi_R(start(j,i):stop(j,i),i);
                    left = pathDB.phi_L(start(j,i):stop(j,i),i);
                    
                    delta = right-left;
                    
                    right = pathDB.phi_R(stop(j,i):start(j+1,i),i);
                    left = pathDB.phi_L(stop(j,i):start(j+1,i),i);
                                        
                    delta2 = right-left;
                    
                    x = linspace(0,1,length(delta)); 
                    x2 = linspace(0,1,length(delta2));
                    
%                     x = frames(start(j,i):stop(j,i));
%                     x2 = frames(stop(j,i):start(j+1,i));
                    
                    y = abs(rad2deg(delta));
                    y2 = abs(rad2deg(delta2));
                    plot(x,y,'color', color(i,1:3))
                    plot(x2,y2, 'color', color(i,1:3))
                end
                end
            end
        end
       
    end
end


ylabel('phi - left minus right [deg]')
hold off

% print(fig, '-djpeg', 'L minus R before stim');