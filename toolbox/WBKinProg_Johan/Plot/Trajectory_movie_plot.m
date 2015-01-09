function Trajectory_movie_plot(settings,pathDB,seq_nr, wingbeat_nr)


    % Program that plots trajectory for Fly_movie4
    
    % Initialize:
    
    % Add path to body model
    
    addpath(char(settings.path_names(7)))
    
    
    addpath(char(settings.path_names(8)))
    
    % Load body model:
    
    cd(char(settings.sequence_names(seq_nr)))
    
    cd('flytracks')
    
    load('ManualFit_flytracks');
    
    cd ..
    
    cd ..
    
    
    % Assign model parameters
    PAR.params = ManualFit.params;
    PAR.DLT = ManualFit.DLT;
    PAR.cam = ManualFit.cam;
    
    clear ManualFit
    
    %--------------------------------
    % Define the Tracking Parameters
    PAR.pixpermm = 1;
    PAR.numfly = 1;
    %Number of parameters of the model (i.e. 8 control points)
    PAR.mdlpar = 15*ones(1,PAR.numfly);
    PAR.statedim = PAR.mdlpar;
    PAR.modelfun_H = @modcurvesplineP;
    PAR.etamax = 0;
    
    %spline order
    PAR.c = 4;
    PAR.L1 = 15; %# of steps for body along length
    PAR.L2 = 6; %# of steps for head along length
    PAR.L3 = 25; %# of steps for wing around the boundary
    PAR.T1 = 13; %# of theta steps for head and body
    PAR.T2 = 2; %# of steps towards center of wing   
    
    
    % Calculate duration wingbeat:
    
    start_meas = find(isnan(pathDB.qb1(:,seq_nr))==0, 1 );
    
    m = wingbeat_nr;
    
       
    
    % Calculate identity begin upstroke:

    
    if pathDB.L_wingbeat_loc(wingbeat_nr,1,seq_nr) > pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)
    
        i_upst_L = pathDB.L_wingbeat_loc(wingbeat_nr,1,seq_nr)-pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)+1;
        
    else
        
        i_upst_L = pathDB.L_wingbeat_loc(wingbeat_nr+1,1,seq_nr)-pathDB.L_wingbeat_loc(wingbeat_nr,2,seq_nr)+1;
        
    end
    
    a = start_meas-1+pathDB.L_wingbeat_loc(m,2,seq_nr)+i_upst_L;
    
    % Calculate current position body:
    
    x_body = pathDB.x_filt(a,seq_nr);
    y_body = pathDB.y_filt(a,seq_nr);
    z_body = pathDB.z_filt(a,seq_nr);
    
    % Calculate trajectory:
    
    x_trajectory = pathDB.x_filt(start_meas:a,seq_nr);
    y_trajectory = pathDB.y_filt(start_meas:a,seq_nr);
    z_trajectory = pathDB.z_filt(start_meas:a,seq_nr);

    % Calculate current body orientation:
    
    q1 = pathDB.qb1_filt(a,seq_nr);
    q2 = pathDB.qb2_filt(a,seq_nr);
    q3 = pathDB.qb3_filt(a,seq_nr);
    q4 = pathDB.qb4_filt(a,seq_nr);
    
    % Create body:
    
    SOLN = [x_body y_body z_body q1 q2 q3 q4 0 0 0 1 0 0 0 1];
    
    clear flymodQ
    [x,y,z] = flymodQ(SOLN,PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    
    % Plot trajectory
    
    hold on
    for k = 1:length(x)
    surf(x{k},y{k},z{k},'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
    end
    plot3(x_trajectory, y_trajectory, z_trajectory,'r')
    axis equal
    hold off
    

end

