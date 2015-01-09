function Trajectory_movie_plot2(settings,pathDB,seq_nr, frame_nr)


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

    start_meas = find(isnan(pathDB.qb1(:,seq_nr))==0, 1 );
    
    % Calculate current position body:
    
    x_body = pathDB.x_filt(frame_nr,seq_nr);
    y_body = pathDB.y_filt(frame_nr,seq_nr);
    z_body = pathDB.z_filt(frame_nr,seq_nr);
    
    % Calculate trajectory:
    
    x_trajectory = pathDB.x_filt(start_meas:frame_nr,seq_nr);
    y_trajectory = pathDB.y_filt(start_meas:frame_nr,seq_nr);
    z_trajectory = pathDB.z_filt(start_meas:frame_nr,seq_nr);

    % Calculate current body orientation:
    
    q1 = pathDB.qb1_filt(frame_nr,seq_nr);
    q2 = pathDB.qb2_filt(frame_nr,seq_nr);
    q3 = pathDB.qb3_filt(frame_nr,seq_nr);
    q4 = pathDB.qb4_filt(frame_nr,seq_nr);
    
    qL1 = pathDB.qL1_filt2(frame_nr,seq_nr);
    qL2 = pathDB.qL2_filt2(frame_nr,seq_nr);
    qL3 = pathDB.qL3_filt2(frame_nr,seq_nr);
    qL4 = pathDB.qL4_filt2(frame_nr,seq_nr);
    
    qR1 = pathDB.qR1_filt2(frame_nr,seq_nr);
    qR2 = pathDB.qR2_filt2(frame_nr,seq_nr);
    qR3 = pathDB.qR3_filt2(frame_nr,seq_nr);
    qR4 = pathDB.qR4_filt2(frame_nr,seq_nr);
    
    % Create body:
    
    SOLN = [x_body y_body z_body q1 q2 q3 q4 qL1 qL2 qL3 qL4 qR1 qR2 qR3 qR4];
    
    clear flymodQ
    [x,y,z] = flymodQ(SOLN,PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end

    
    % Plot trajectory
    
    hold on
    for k = 1:length(x)
    surf(x{k},y{k},z{k},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong');
    end
    plot3(x_trajectory, y_trajectory, z_trajectory,'r')
    axis equal
    hold off
    

end