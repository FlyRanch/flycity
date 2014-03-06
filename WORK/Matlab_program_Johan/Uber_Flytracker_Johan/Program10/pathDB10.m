function pathDB10( settings, pathDB )

    savefile = 'pathDB10';
    
    
    [control_gen,control_glob] = control_analysis( pathDB );
    
    save(savefile,'control_glob','control_gen')


end

