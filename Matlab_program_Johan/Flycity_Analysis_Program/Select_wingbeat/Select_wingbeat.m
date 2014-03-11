function Select_wingbeat(settings,pathDB)

    savefile = 'pathDB4.mat';

    % Select the beginning and the end of the wingbeat
    
    N = settings.nr_of_seq;
    
    M = settings.frame_end;
    
    wingbeats.nr_of_wb              = nan(N,1);
    wingbeats.wingbeat_duration     = nan(200,N);
    wingbeats.downstroke_duration   = nan(200,N);
    wingbeats.upstroke_duration     = nan(200,N);
    wingbeats.wingbeat_loc          = nan(200,2,N);
    wingbeats.downstroke_loc        = nan(200,2,N);
    wingbeats.upstroke_loc          = nan(200,2,N);
    
    t = pathDB.t;
    
    % Construct the wingtip path:
    
    wingtip_path_L = nan(3,M,N);
    wingtip_path_R = nan(3,M,N);
    
    for i = 1:N
               
        for j = 1:M
            
            RL = pathDB.rot_mat.RL(:,:,j,i);
            RR = pathDB.rot_mat.RR(:,:,j,i);
            
            wingtip_path_L(:,j,i) = RL*[0; -1; 0];
            wingtip_path_R(:,j,i) = RR*[0; 1; 0];
            
        end
        
    end
    
    
    % Find the upstroke and downstroke location:
    
    for i = 1:N
        
        % Find the locations of the start of the upstroke and the start of the
        % downstroke (sequences will start with the first downstroke)
        
        start = settings.start_stop(i,1);
        stop = settings.start_stop(i,2);
        
        [ GammaL, locL] = findpeaks(wingtip_path_L(2,:,i), 'minpeakdistance' , 8);
        [ GammaR, locR] = findpeaks((-wingtip_path_R(2,:,i)), 'minpeakdistance' , 8);
        
        loc_match = nan(200,1);
                        
        for j = 1:length(locL)
            
            if locR(j) <= locL(j)+3 && locR(j) >= locL(j)-3
                
                loc_match(j) = 1;
                
            else
                
                loc_match(j) = 0;
                
            end
            
        end
        
        nr_wb = find(isnan(loc_match)==1,1,'first')-1;
        
        loc_match_zeros = find(loc_match == 0);
        
        if isempty(loc_match_zeros) == 0
                
            if loc_match_zeros(1) == 1
                
                if locL(1) <= start+2 && loc_match(2) == 1
                    
                    loc_match(1) = 1;
                    loc_match_zeros(1) = [];
                    
                end
                
            end
                
            if loc_match_zeros(end) == nr_wb
                
                if locL(nr_wb) >= stop-2 && loc_match(nr_wb-1) == 1
                    
                    loc_match(nr_wb) = 1;
                    loc_match_zeros(nr_wb) = [];
                    
                end
                
            end
                
            if isempty(loc_match_zeros) == 0
                
                if length(loc_match_zeros) == 1
                    
                    wb_intervals = [ loc_match_zeros(1)-1 nr_wb-loc_match_zeros(1) ];
                    
                else
                    
                    wb_intervals = [ loc_match_zeros(1)-1 loc_match_zeros(2:end)-loc_match_zeros(1:(end-1))-1 nr_wb-loc_match_zeros(end) ];
                    
                end
                
                [max_interval max_loc ] = max(wb_intervals);
                
                if max_loc == 1
                    
                    loc_match((loc_match_zeros(1)):end) = nan;
                    
                elseif max_loc == length(wb_intervals)
                                        
                    loc_match(1:loc_match_zeros(max_loc-1)) = nan;
                    
                else
                    
                    loc_match(1:(loc_match_zeros(max_loc-1))) = nan;
                    loc_match(loc_match_zeros(max_loc):end) = nan;
                    
                end
                
            end
            
                            
        end
        
        wb_start = find(loc_match==1, 1, 'first');
        wb_end   = find(isnan(loc_match)==1, 1, 'first')-1;
        
        if abs(sum(GammaL(wb_start:2:wb_end))) > abs(sum(GammaL((wb_start+1):2:wb_end)))

            loc_upst_L = locL(wb_start:2:wb_end);
            loc_dwnst_L = locL((wb_start+1):2:wb_end);

        elseif abs(sum(GammaL(wb_start:2:wb_end))) < abs(sum(GammaL((wb_start+1):2:wb_end)))

            loc_upst_L = locL((wb_start+1):2:wb_end);
            loc_dwnst_L = locL(wb_start:2:wb_end);

        end
        
                
        if loc_upst_L(1) < loc_dwnst_L(1)
            
            loc_upst_L(1) = [];
            
        end
        
        if loc_upst_L(end) < loc_dwnst_L(end)
            
            last_wb_end = loc_dwnst_L(end)-1;
            
            loc_dwnst_L(end) = [];
            
        elseif loc_upst_L(end) > loc_dwnst_L(end)
            
            loc_upst_L(end) = [];
            
            last_wb_end = loc_dwnst_L(end)-1;
            
            loc_dwnst_L(end) = [];
                       
        end
        
        if length(loc_dwnst_L) ~= length(loc_upst_L)
            
            'number of up and down strokes is not equal'
            
        end
        
        nr_wb2 = length(loc_dwnst_L);
        
        wingbeats.nr_of_wb(i)                     = nr_wb2;
        wingbeats.wingbeat_loc(1:nr_wb2,:,i)      = [ loc_dwnst_L(1:nr_wb2)' [(loc_dwnst_L(2:nr_wb2)-1)'; last_wb_end] ];
        wingbeats.downstroke_loc(1:nr_wb2,:,i)   = [ loc_dwnst_L(1:nr_wb2)' (loc_upst_L(1:nr_wb2)-1)' ];
        wingbeats.upstroke_loc(1:nr_wb2,:,i)     = [ loc_upst_L(1:nr_wb2)' [(loc_dwnst_L(2:nr_wb2)-1)'; last_wb_end] ];
        wingbeats.wingbeat_duration(1:nr_wb2,i)   = wingbeats.wingbeat_loc(1:nr_wb2,2,i)-wingbeats.wingbeat_loc(1:nr_wb2,1,i)+1;
        wingbeats.downstroke_duration(1:nr_wb2,i) = wingbeats.downstroke_loc(1:nr_wb2,2,i)-wingbeats.downstroke_loc(1:nr_wb2,1,i)+1;
        wingbeats.upstroke_duration(1:nr_wb2,i)   = wingbeats.upstroke_loc(1:nr_wb2,2,i)-wingbeats.upstroke_loc(1:nr_wb2,1,i)+1;
        
%         figure()
%         hold on
%         plot(wingtip_path_L(2,:,i))
%         plot(locL,wingtip_path_L(2,locL,i),'o','Color','r')
%         hold off
%         
%         figure()
%         hold on
%         plot(-wingtip_path_R(2,:,i))
%         plot(locR,-wingtip_path_R(2,locR,i),'o','Color','r')
%         hold off
%         
%         figure(i)
%         hFig = figure(i);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 1400 800]);
%         hold on
%         plot(t(start:stop),wingtip_path_L(2,start:stop,i))
%         plot(t(loc_dwnst_L),wingtip_path_L(2,loc_dwnst_L,i),'o','Color','r')
%         plot(t(loc_upst_L),wingtip_path_L(2,loc_upst_L,i),'o','Color','g')
%         hold off
%         title('Distance left wingtip - symmetry plane')
%         xlabel('t [s]')
%         ylabel('y_L [mm]')
%         
%         pause
        
    end
    
    
   
    % Save the structures:
    
    save(savefile,'wingbeats')

end

