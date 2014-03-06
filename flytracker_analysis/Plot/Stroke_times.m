function [ upstroke_t, downstroke_t ] = Stroke_times( wingtip_y, wingbeat_loc, time)

    % Program that calculates the begin and end time of the upstroke and
    % the downstroke and returns two vectors with the begin and end time:
    
    % start and stop point for the measurements
    start = find(isnan(wingbeat_loc(:))==0, 1 );
    stop = find(isnan(wingbeat_loc(:))==0, 1, 'last' );
    
    N = stop-start+1;
    
    upstroke_t = zeros(2*ceil(N/2),1);
    downstroke_t = zeros(2*ceil(N/2),1);
    
    j = 1;
    
    for i = 1:N-1
        
        if wingtip_y(wingbeat_loc(i))>wingtip_y(wingbeat_loc(i+1))
            
            downstroke_t(j) = time(wingbeat_loc(i));
            downstroke_t(j+1) = time(wingbeat_loc(i+1));
            
        elseif wingtip_y(wingbeat_loc(i))<wingtip_y(wingbeat_loc(i+1))
            
            upstroke_t(j) = time(wingbeat_loc(i));
            upstroke_t(j+1) = time(wingbeat_loc(i+1));
            
        else
            
        %    ' upstroke-downstroke determination fails'
            
        end
        
        j = j+2;
        
    end
    
    
end