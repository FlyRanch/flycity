function Fly_movie_forces( body_pos, Rb, RL, RR, F_tot, FA, FI, Fg, body_model, wing_model, movie_name,  )


    % Initialize movie:
    
    movfilename = [ 'force_movie' char(movie_name) '.mj2' ];
    mov = VideoWriter(movfilename,'Motion JPEG AVI');
    mov.FrameRate = 30;
    mov.Quality = 100;
    
    N = size(RL,3);
    
    for i = 1:N
        
        
    
        Fly_plot_3D( body_pos, Rb, RL, RR, body_model, wing_model )
    
    end

end