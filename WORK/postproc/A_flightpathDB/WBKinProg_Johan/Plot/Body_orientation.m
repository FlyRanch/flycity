function Body_orientation(settings, pathDB, seq_nr, phi, theta, xsi, frame_nr)


    % Program that generates 3D plots of the wingtip_path, postion of the
    % fruit fly body and the orientation of the wing for a single given
    % wingbeat
    
    % Initialize:
    
    
    % Load body model:
    
    cd(char(settings.sequence_names(seq_nr)))
    
    addpath(char(settings.path_names(7)))
    
    addpath(char(settings.path_names(8)))
    
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


  
    R_phi   = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R_theta = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R_xsi   = [cos(xsi) -sin(xsi) 0; sin(xsi) cos(xsi) 0; 0 0 1];
    
    R_matrix = R_phi*R_theta*R_xsi;
    
    q_body = [sin(phi/2)*cos(theta/2)*cos(xsi/2)-cos(phi/2)*sin(theta/2)*sin(xsi/2); ...
              cos(phi/2)*sin(theta/2)*cos(xsi/2)+sin(phi/2)*cos(theta/2)*sin(xsi/2); ...
              cos(phi/2)*cos(theta/2)*sin(xsi/2)-sin(phi/2)*sin(theta/2)*cos(xsi/2); ...
              cos(phi/2)*cos(theta/2)*cos(xsi/2)+sin(phi/2)*sin(theta/2)*sin(xsi/2)];
          
    q_body = q_body./norm(q_body);
    
    q_wing_L = [pathDB.qL1_filt2(frame_nr,seq_nr); pathDB.qL2_filt2(frame_nr,seq_nr); pathDB.qL3_filt2(frame_nr,seq_nr); pathDB.qL4_filt2(frame_nr,seq_nr)];
    
    q_wing_R = [pathDB.qR1_filt2(frame_nr,seq_nr); pathDB.qR2_filt2(frame_nr,seq_nr); pathDB.qR3_filt2(frame_nr,seq_nr); pathDB.qR4_filt2(frame_nr,seq_nr)];
    
    SOLN = [0 0 0 q_body(1) q_body(2) q_body(3) q_body(4) q_wing_L(1) q_wing_L(2) q_wing_L(3) q_wing_L(4) q_wing_R(1) q_wing_R(2) q_wing_R(3) q_wing_R(4)];
    
    clear flymodQ
    [x,y,z] = flymodQ(SOLN,PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end
    
    hold on
    for k = 1:length(x)
    surf(x{k},y{k},z{k},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong');
    end
    axis equal
    hold off



end

