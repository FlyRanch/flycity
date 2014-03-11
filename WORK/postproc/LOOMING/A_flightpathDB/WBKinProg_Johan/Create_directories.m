function [plot_folders] = Create_directories( settings )

    % Program that creates the folders in which plots and data of the film
    % sequences will be stored.
    
    % Create a string including date and time
    
    date_time = clock;
    
    time_str = [ 'Plots_' int2str(date_time(3)) '_' int2str(date_time(2)) '_' int2str(date_time(1)) '_' int2str(date_time(4)) '_' int2str(date_time(5)) '_' int2str(date_time(6)) ];
    
    % Create a folder in the folder Plots which will contain all the plots
    % and data of the current run:
    
    cd(char(settings.path_names(6)))
    
    mkdir(time_str)
    
    cd([char(settings.path_names(6)) '/' time_str])
    
    % Create the folder sequences
    
    mkdir('sequences')
    
    % Create all the sequences in this folder:
    
    cd([char(settings.path_names(6)) '/' time_str '/sequences'])
    
    for j = 1:length(settings.sequence_names)
        
        mkdir(char(settings.sequence_names(j)))
        
    end
    
    cd([char(settings.path_names(6)) '/' time_str])
    
    % Create the folder trajectories:
    
    mkdir('trajectories')
    
    traject_names = {'flight_path'; ...
                     'xyz_body'; ...
                     'uvw_body'; ...
                     'axayaz_body'; ...
                     'quat_body'; ...
                     'w_body'; ...
                     'ypr_body'; ...
                     'quat_left'; ...
                     'w_left'; ...
                     'pte_left'; ...
                     'quat_right'; ...
                     'w_right'; ...
                     'pte_right'; ...
                     'pte_bref_lr'; ...
                     'sphere_plot'; ...
                     'pte_stroke_lr'; ...
                     'alfa_beta_body'; ...
                     'body_vel'; ...
                     'wing_vel'; ...
                     'alfa_beta_wing'; ...
                     'trans_rot_wing'};
    
    wingbeat_names = {'Wing_kinematics_3D_stationary'; ...
                      'Wing_kinematics_3D_moving'; ...
                      'Wing_kinematics_2D'; ...
                      'Arm_force'; ...
                      'Force_3D_stationary'; ...
                      'Force_3D_moving'};
                 
    % Create all trajectories in this folder:
    
    cd([char(settings.path_names(6)) '/' time_str '/trajectories'])
                 
    for j = 1:length(traject_names)
        
        mkdir(char(traject_names(j)))
        
    end
    
    cd([char(settings.path_names(6)) '/' time_str])
    
    % Create the folder wingbeats:
           
    mkdir('wingbeats')
    
    cd([char(settings.path_names(6)) '/' time_str '/wingbeats'])
    
    for j = 1:length(settings.sequence_names)
        
        mkdir(char(settings.sequence_names(j)))
        
        cd(char(settings.sequence_names(j)))
        
        for k = 1:length(wingbeat_names)
        
        mkdir(char(wingbeat_names(k)))
        
        end
        
        cd ..
    end
    
    cd([char(settings.path_names(6)) '/' time_str])
    
    % Create the folder movies:
    
    mkdir('movies')
    
    cd([char(settings.path_names(6)) '/' time_str '/movies'])
    
    for j = 1:length(settings.sequence_names)
        
        mkdir(char(settings.sequence_names(j)))
        
    end
    
    % Return the folder locations:
    
    plot_folders = {[ char(settings.path_names(6)) '/' time_str '/sequences']; ...
                    [ char(settings.path_names(6)) '/' time_str '/trajectories']; ...
                    [ char(settings.path_names(6)) '/' time_str '/wingbeats']; ...
                    [ char(settings.path_names(6)) '/' time_str '/movies']};
                
    % Reset current directory to original value:
    
    cd(char(settings.path_names(1)))
                    

end

