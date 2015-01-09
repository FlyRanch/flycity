function [ x, y, z ] = load_body_model(settings, seq_nr, xh )

    
    addpath(char(settings.path_names(7)))
    addpath(char(settings.path_names(8)))

    % Program that loads and creates body model for certain set of
    % quaternions

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
    
    clear flymodQ
    [x,y,z] = flymodQ(xh,PAR.params,PAR);
    for j = 1:length(x);
        PAR.modsample(j) = size(x{j},1);
    end


end

