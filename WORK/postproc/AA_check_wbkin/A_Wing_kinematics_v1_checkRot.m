%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Created at: October, 19th, 2012
%                   By: Johan Melis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all
% clc

if exist('WBkin.mat') == 2

    load('WBkin.mat')

else

%% constants
% strokeplane_WBkin = -55;
% strokeplane_WBkin = -45;
strokeplane_WBkin = -47.5;

settings.strokeplane_WBkin = strokeplane_WBkin;

%%

add_paths = {...
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan'; ...
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan/Kalman_Filter'; ...
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan/Plot'; ...
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan/Plot/Smoothing_plots'; ... % };
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan/Body_model'; ...
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan/Body_model/mex'; ...
             '/home/florian/Dropbox/FlyCity/WORK/postproc/AA_check_wbkin/WBKinProg_Johan/Type2Regression';...
             };
         
addpath(char(add_paths(1)))
addpath(char(add_paths(2)))
addpath(char(add_paths(3)))
addpath(char(add_paths(4)))
         
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/');
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/mex/');
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/core/');
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/results/');

         
% add_paths = {'F:/Sequences_ABCD/A'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Kalman_Filter'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot/Smoothing_plots'; ... % };
%              'F:/Sequences_ABCD/Plots'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model/mex'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Type2Regression';...
%              };
% 
% add_paths = {'F:/wingsInc/All_tracked_seq'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Kalman_Filter'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot/Smoothing_plots'; ... % };
%              'F:/wingsInc/Plots'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model/mex'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Type2Regression';...
%              };

% add_paths = {'D:/FLIGHT_SEQS/complete'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Kalman_Filter'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot/Smoothing_plots'; ... % };
%              'D:/FLIGHT_SEQS/Plots'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model/mex'; ...
%              'D:/FLIGHT_SEQS/Projections_20121010_1'};

% add_paths = {'F:/FLIGHT_SEQS/complete'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Kalman_Filter'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Plot/Smoothing_plots'; ... % };
%              'F:/FLIGHT_SEQS/Plots'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model'; ...
%              'D:/Flytracker/Wing_kin_Johan/Program9/Body_model/mex'; ...
%              'F:/FLIGHT_SEQS/Projections_20121010_1'};
%              
%             %              'D:/FlyData/Plots'};
% 
% 
% addpath(char(add_paths(1)))
% addpath(char(add_paths(2)))
% addpath(char(add_paths(3)))
% addpath(char(add_paths(4)))




%%

% initialize settings
settings.path_names = add_paths;
settings.fps = 7500;
settings.frame_end = 5588;
settings.trigger_frame = 2795;


% Determine whether data is going to be saved or not [0 = off, 1 = on]:
settings.save_on_off = 0;

temp_dir = dir;
root_dir = cd;

% Add folder names containing flight data to settings.folder_names

% if length(temp_dir)<3
%     
%     'no data files in directory'
%     pause
%     
% else
%     
%     
%     temp_dirnames = {temp_dir.name};
%     
%     % FTMmod: change folder structure
%     temp_foldernames = []
%     for i = 3:length(temp_dir)
%          dir_now = temp_dir(i);
%        
%         if dir_now.isdir == 1
%             cd(dir_now.name)
%             seq_dir = dir('*S00*');
%             
%             if ~isempty(seq_dir)
%                 for j = 1:length(seq_dir)
% 
%                     temp_seqname = [cd,'/',seq_dir(j).name];
%                     temp_foldernames = [temp_foldernames; {temp_seqname}];
%                 end
%             end
%             cd ..
%         end
%     end
% end
    settings.sequence_names = [{cd}];
%     settings.root_dir = root_dir;
            
        
%         b = cell2mat(strfind(temp_dirnames(i),'S00'));
%         
%         if b == 10
%             
%             temp_foldernames = [temp_foldernames; temp_dirnames(i)];
%             
%         end
%         clear b
%     end
%     
%     settings.sequence_names = temp_foldernames;
% 
% end
% 
if settings.save_on_off == 1
    settings.plot_folders = Create_directories( settings );
end


clear temp_dirnames temp_dir temp_foldernames i



% initialize DB variable

% Computation variables

frames = [1:settings.frame_end]';
pathDB.t = (frames-settings.trigger_frame)/settings.fps;

pathDB.x = [];
pathDB.y = [];
pathDB.z = [];

pathDB.x_filt = [];
pathDB.y_filt = [];
pathDB.z_filt = [];

pathDB.u_filt = [];
pathDB.v_filt = [];
pathDB.w_filt = [];

pathDB.ax_filt = [];
pathDB.ay_filt = [];
pathDB.az_filt = [];

pathDB.u_body = [];
pathDB.v_body = [];
pathDB.w_body = [];

pathDB.ax_body = [];
pathDB.ay_body = [];
pathDB.az_body = [];

pathDB.qb1 = [];
pathDB.qb2 = [];
pathDB.qb3 = [];
pathDB.qb4 = [];

pathDB.qL1 = [];
pathDB.qL2 = [];
pathDB.qL3 = [];
pathDB.qL4 = [];

pathDB.qR1 = [];
pathDB.qR2 = [];
pathDB.qR3 = [];
pathDB.qR4 = [];

pathDB.qb1_filt = [];
pathDB.qb2_filt = [];
pathDB.qb3_filt = [];
pathDB.qb4_filt = [];

pathDB.b_omega1 = [];
pathDB.b_omega2 = [];
pathDB.b_omega3 = [];

pathDB.qL1_filt1 = [];
pathDB.qL2_filt1 = [];
pathDB.qL3_filt1 = [];
pathDB.qL4_filt1 = [];

pathDB.qR1_filt1 = [];
pathDB.qR2_filt1 = [];
pathDB.qR3_filt1 = [];
pathDB.qR4_filt1 = [];

pathDB.omega1_L = [];
pathDB.omega2_L = [];
pathDB.omega3_L = [];

pathDB.qL1_filt2 = [];
pathDB.qL2_filt2 = [];
pathDB.qL3_filt2 = [];
pathDB.qL4_filt2 = [];

pathDB.omega1_R = [];
pathDB.omega2_R = [];
pathDB.omega3_R = [];

pathDB.qR1_filt2 = [];
pathDB.qR2_filt2 = [];
pathDB.qR3_filt2 = [];
pathDB.qR4_filt2 = [];

pathDB.b_roll = [];
pathDB.b_pitch = [];
pathDB.b_yaw = [];
pathDB.b_alfa = [];
pathDB.b_beta = [];

pathDB.body_l = [];
pathDB.wing_l = [];
pathDB.joint_pos_L = [];
pathDB.joint_pos_R = [];

pathDB.u_wing_L = [];
pathDB.v_wing_L = [];
pathDB.w_wing_L = [];

pathDB.u_wing_R = [];
pathDB.v_wing_R = [];
pathDB.w_wing_R = [];

pathDB.u_wing_rot_L = [];
pathDB.v_wing_rot_L = [];
pathDB.w_wing_rot_L = [];
    
pathDB.u_wing_rot_R = [];
pathDB.v_wing_rot_R = [];
pathDB.w_wing_rot_R = [];
    
pathDB.u_wing_trans_L = [];
pathDB.v_wing_trans_L = [];
pathDB.w_wing_trans_L = [];
    
pathDB.u_wing_trans_R = [];
pathDB.v_wing_trans_R = [];
pathDB.w_wing_trans_R = [];
    
pathDB.u_wing_body_rot_L = [];
pathDB.v_wing_body_rot_L = [];
pathDB.w_wing_body_rot_L = [];
    
pathDB.u_wing_body_rot_R = [];
pathDB.v_wing_body_rot_R = [];
pathDB.w_wing_body_rot_R = [];
    
pathDB.alfa_L = [];
pathDB.beta_L = [];
    
pathDB.alfa_R = [];
pathDB.beta_R = [];

pathDB.Lwingtip = [];
pathDB.Rwingtip = [];

pathDB.L_wingbeat_loc = [];
pathDB.R_wingbeat_loc = [];

pathDB.phi_L = [];
pathDB.theta_L = [];
pathDB.eta_L = [];

pathDB.phi_R = [];
pathDB.theta_R = [];
pathDB.eta_R = [];

% Histogram data:

    % Body:
    
    pathDB.phi_body_mean = [];
    pathDB.theta_body_mean = [];
    pathDB.xsi_body_mean = [];
    
    pathDB.phi_body_sd = [];
    pathDB.theta_body_sd = [];
    pathDB.xsi_body_sd = [];
    
    pathDB.omegax_body_mean = [];
    pathDB.omegay_body_mean = [];
    pathDB.omegaz_body_mean = [];
    pathDB.Omega_body_mean = [];
    
    pathDB.omegax_body_sd = [];
    pathDB.omegay_body_sd = [];
    pathDB.omegaz_body_sd = [];
    pathDB.Omega_body_sd = [];
    
    pathDB.alfa_body_mean = [];
    pathDB.beta_body_mean = [];
    
    pathDB.alfa_body_sd = [];
    pathDB.beta_body_sd = [];
    
    pathDB.u_body_mean = [];
    pathDB.v_body_mean = [];
    pathDB.w_body_mean = [];
    pathDB.U_body_mean = [];
    
    pathDB.u_body_sd = []; 
    pathDB.v_body_sd = [];
    pathDB.w_body_sd = [];
    pathDB.U_body_sd = [];
    
    pathDB.ax_body_mean = [];
    pathDB.ay_body_mean = [];
    pathDB.az_body_mean = [];
    pathDB.a_body_mean = [];
    
    pathDB.ax_body_sd = [];
    pathDB.ay_body_sd = [];
    pathDB.az_body_sd = [];
    pathDB.a_body_sd = [];
    
    % Left wing:
    
    pathDB.phi_L_up_mean = []; 
    pathDB.theta_L_up_mean = [];
    pathDB.eta_L_up_mean = [];
    
    pathDB.phi_L_up_sd = [];
    pathDB.theta_L_up_sd = [];
    pathDB.eta_L_up_sd = [];
    
    pathDB.phi_L_down_mean = [];
    pathDB.theta_L_down_mean = [];
    pathDB.eta_L_down_mean = [];
    
    pathDB.phi_L_down_sd = [];
    pathDB.theta_L_down_sd = [];
    pathDB.eta_L_down_sd = [];
    
    pathDB.alfa_L_up_mean = []; 
    pathDB.beta_L_up_mean = [];
    
    pathDB.alfa_L_up_sd = [];
    pathDB.beta_L_up_sd = [];
    
    pathDB.alfa_L_down_mean = [];
    pathDB.beta_L_down_mean = [];
    
    pathDB.alfa_L_down_sd = [];
    pathDB.beta_L_down_sd = [];
    
    pathDB.u_L_up_mean = [];
    pathDB.v_L_up_mean = [];
    pathDB.w_L_up_mean = [];
    pathDB.U_L_up_mean = [];
    
    pathDB.u_L_up_sd = [];
    pathDB.v_L_up_sd = [];
    pathDB.w_L_up_sd = [];
    pathDB.U_L_up_sd = [];
    
    pathDB.u_L_down_mean = [];
    pathDB.v_L_down_mean = [];
    pathDB.w_L_down_mean = [];
    pathDB.U_L_down_mean = [];
    
    pathDB.u_L_down_sd = []; 
    pathDB.v_L_down_sd = [];
    pathDB.w_L_down_sd = [];
    pathDB.U_L_down_sd = [];
    
    pathDB.omegax_L_up_mean = [];
    pathDB.omegay_L_up_mean = [];
    pathDB.omegaz_L_up_mean = [];
    pathDB.Omega_L_up_mean = [];
    
    pathDB.omegax_L_up_sd = [];
    pathDB.omegay_L_up_sd = [];
    pathDB.omegaz_L_up_sd = [];
    pathDB.Omega_L_up_sd = [];
    
    pathDB.omegax_L_down_mean = [];
    pathDB.omegay_L_down_mean = [];
    pathDB.omegaz_L_down_mean = [];
    pathDB.Omega_L_down_mean = [];
    
    pathDB.omegax_L_down_sd = [];
    pathDB.omegay_L_down_sd = [];
    pathDB.omegaz_L_down_sd = [];
    pathDB.Omega_L_down_sd = [];
    
    
    % Right wing:
    
    pathDB.phi_R_up_mean = [];
    pathDB.theta_R_up_mean = [];
    pathDB.eta_R_up_mean = [];
    
    pathDB.phi_R_up_sd = [];
    pathDB.theta_R_up_sd = [];
    pathDB.eta_R_up_sd = [];
    
    pathDB.phi_R_down_mean = [];
    pathDB.theta_R_down_mean = [];
    pathDB.eta_R_down_mean = [];
    
    pathDB.phi_R_down_sd = [];
    pathDB.theta_R_down_sd = [];
    pathDB.eta_R_down_sd = [];
    
    pathDB.alfa_R_up_mean = [];
    pathDB.beta_R_up_mean = [];
    
    pathDB.alfa_R_up_sd = [];
    pathDB.beta_R_up_sd = [];
    
    pathDB.alfa_R_down_mean = [];
    pathDB.beta_R_down_mean = [];
    
    pathDB.alfa_R_down_sd = [];
    pathDB.beta_R_down_sd = [];
    
    pathDB.u_R_up_mean = [];
    pathDB.v_R_up_mean = [];
    pathDB.w_R_up_mean = [];
    pathDB.U_R_up_mean = [];
    
    pathDB.u_R_up_sd = [];
    pathDB.v_R_up_sd = [];
    pathDB.w_R_up_sd = [];
    pathDB.U_R_up_sd = [];
    
    pathDB.u_R_down_mean = [];
    pathDB.v_R_down_mean = [];
    pathDB.w_R_down_mean = [];
    pathDB.U_R_down_mean = [];
    
    pathDB.u_R_down_sd = [];
    pathDB.v_R_down_sd = [];
    pathDB.w_R_down_sd = [];
    pathDB.U_R_down_sd = [];
    
    pathDB.omegax_R_up_mean = [];
    pathDB.omegay_R_up_mean = [];
    pathDB.omegaz_R_up_mean = [];
    pathDB.Omega_R_up_mean = [];
    
    pathDB.omegax_R_up_sd = [];
    pathDB.omegay_R_up_sd = [];
    pathDB.omegaz_R_up_sd = [];
    pathDB.Omega_R_up_sd = [];
    
    pathDB.omegax_R_down_mean = [];
    pathDB.omegay_R_down_mean = [];
    pathDB.omegaz_R_down_mean = [];
    pathDB.Omega_R_down_mean = [];
    
    pathDB.omegax_R_down_sd = [];
    pathDB.omegay_R_down_sd = [];
    pathDB.omegaz_R_down_sd = [];
    pathDB.Omega_R_down_sd = [];

    
    % PathDB7
   
    pathDB.F_joint_L_down = [];
    pathDB.M_joint_L_down = [];
    
    pathDB.F_joint_L_up = [];
    pathDB.M_joint_L_up = [];
    
    pathDB.F_joint_R_down = [];
    pathDB.M_joint_R_down = [];
    
    pathDB.F_joint_R_up = [];
    pathDB.M_joint_R_up = [];
    
    pathDB.F_cg = [];
    pathDB.M_cg = [];
    
    pathDB.R_L_down = [];
    pathDB.Fn_L_down = [];
    pathDB.Ft_L_down = [];
    pathDB.Mn_L_down = [];
    pathDB.Mt_L_down = [];
    
    pathDB.R_L_up = [];
    pathDB.Fn_L_up = [];
    pathDB.Ft_L_up = [];
    pathDB.Mn_L_up = [];
    pathDB.Mt_L_up = [];
    
    pathDB.R_R_down = [];
    pathDB.Fn_R_down = [];
    pathDB.Ft_R_down = [];
    pathDB.Mn_R_down = [];
    pathDB.Mt_R_down = [];
    
    pathDB.R_R_up = [];
    pathDB.Fn_R_up = [];
    pathDB.Ft_R_up = [];
    pathDB.Mn_R_up = [];
    pathDB.Mt_R_up = [];
      
    pathDB.m_est_point = [];
    
    pathDB.m_est_dynamical_model = [];
    
    pathDB.r_body_L_down = [];
    pathDB.r_body_L_up = [];
    
    pathDB.r_body_R_down = [];
    pathDB.r_body_R_up = [];
    
    pathDB.q_avg_down_body = [];
    pathDB.q_avg_up_body = [];
    
    pathDB.q_avg_down_L = [];
    pathDB.q_avg_up_L = [];
    
    pathDB.q_avg_down_R = [];
    pathDB.q_avg_up_R = [];
    
    pathDB.down_time_L = [];
    pathDB.up_time_L = [];
        
    pathDB.down_time_R= [];
    pathDB.up_time_R= [];
    
    pathDB.wingbeat_time= [];
    
    pathDB.F_L = [];
    pathDB.F_R = [];
    
    pathDB.F_L_arm = [];
    pathDB.F_R_arm = [];
    
    
% PathDB8


    

%%

% load or create a file with the observed x,y,z data

if exist('pathDB1.mat')==2
    temp = load('pathDB1.mat');
        pathDB.x = temp.x;
        pathDB.y = temp.y;
        pathDB.z = temp.z;

        pathDB.qb1 = temp.bodyQ1;
        pathDB.qb2 = temp.bodyQ2;
        pathDB.qb3 = temp.bodyQ3;
        pathDB.qb4 = temp.bodyQ4;

        pathDB.qL1 = temp.wingQL1;
        pathDB.qL2 = temp.wingQL2;
        pathDB.qL3 = temp.wingQL3;
        pathDB.qL4 = temp.wingQL4;

        pathDB.qR1 = temp.wingQR1;
        pathDB.qR2 = temp.wingQR2;
        pathDB.qR3 = temp.wingQR3;
        pathDB.qR4 = temp.wingQR4;
        
    'loaded'
    clear temp
else
    temp = make_pathDB1_nosave(settings,pathDB);
%     temp = load('pathDB1.mat');
        pathDB.x = temp.x;
        pathDB.y = temp.y;
        pathDB.z = temp.z;

        pathDB.qb1 = temp.bodyQ1;
        pathDB.qb2 = temp.bodyQ2;
        pathDB.qb3 = temp.bodyQ3;
        pathDB.qb4 = temp.bodyQ4;

        pathDB.qL1 = temp.wingQL1;
        pathDB.qL2 = temp.wingQL2;
        pathDB.qL3 = temp.wingQL3;
        pathDB.qL4 = temp.wingQL4;

        pathDB.qR1 = temp.wingQR1;
        pathDB.qR2 = temp.wingQR2;
        pathDB.qR3 = temp.wingQR3;
        pathDB.qR4 = temp.wingQR4;
        
%     'created + loaded'
    clear temp
end


%%
% 
% %Filter the body position and body orientation with a Kalman filter.
% 
% if exist('pathDB2.mat')==2
%     temp = load('pathDB2.mat');
%         pathDB.x_filt = temp.x;
%         pathDB.y_filt = temp.y;
%         pathDB.z_filt = temp.z;
%         
%         pathDB.u_filt = temp.u;
%         pathDB.v_filt = temp.v;
%         pathDB.w_filt = temp.w;
% 
%         pathDB.ax_filt = temp.ax;
%         pathDB.ay_filt = temp.ay;
%         pathDB.az_filt = temp.az;
% 
%         pathDB.qb1_filt = temp.q1;
%         pathDB.qb2_filt = temp.q2;
%         pathDB.qb3_filt = temp.q3;
%         pathDB.qb4_filt = temp.q4;
%         
%         pathDB.b_omega1 = temp.omega_1;
%         pathDB.b_omega2 = temp.omega_2;
%         pathDB.b_omega3 = temp.omega_3;
%         
%     'loaded'
%     clear temp
% else
%     temp = filter_body_nosave(settings,pathDB);
% %     temp = load('pathDB2.mat');
%         pathDB.x_filt = temp.x;
%         pathDB.y_filt = temp.y;
%         pathDB.z_filt = temp.z;
%         
%         pathDB.u_filt = temp.u;
%         pathDB.v_filt = temp.v;
%         pathDB.w_filt = temp.w;
% 
%         pathDB.ax_filt = temp.ax;
%         pathDB.ay_filt = temp.ay;
%         pathDB.az_filt = temp.az;
% 
%         pathDB.qb1_filt = temp.q1;
%         pathDB.qb2_filt = temp.q2;
%         pathDB.qb3_filt = temp.q3;
%         pathDB.qb4_filt = temp.q4;
%         
%         pathDB.b_omega1 = temp.omega_1;
%         pathDB.b_omega2 = temp.omega_2;
%         pathDB.b_omega3 = temp.omega_3;
%         
%         
% %     'created + loaded'
%     clear temp
% end
% 

%%

% Correct wing quaternions for filtered body quaternion and calculate the
% body angles psi, theta, phi, alfa and beta and their derivatives.

if exist('pathDB3.mat')==2
    temp = load('pathDB3.mat');
    
    pathDB.qL1_filt1 = temp.qL1_filt1;
    pathDB.qL2_filt1 = temp.qL2_filt1;
    pathDB.qL3_filt1 = temp.qL3_filt1;
    pathDB.qL4_filt1 = temp.qL4_filt1;

    pathDB.qR1_filt1 = temp.qR1_filt1;
    pathDB.qR2_filt1 = temp.qR2_filt1;
    pathDB.qR3_filt1 = temp.qR3_filt1;
    pathDB.qR4_filt1 = temp.qR4_filt1;
    
    pathDB.b_alfa = temp.alfa_body;
    pathDB.b_beta = temp.beta_body;
    
    pathDB.b_roll = temp.roll_body;
    pathDB.b_pitch = temp.pitch_body;
    pathDB.b_yaw = temp.yaw_body;
    
    pathDB.u_body = temp.u_body;
    pathDB.v_body = temp.v_body;
    pathDB.w_body = temp.w_body;
    
    pathDB.ax_body = temp.ax_body;
    pathDB.ay_body = temp.ay_body;
    pathDB.az_body = temp.az_body;
    
    pathDB.body_l = temp.body_length;
    pathDB.wing_l = temp.wing_length;
    pathDB.joint_pos_L = temp.joint_pos_L;
    pathDB.joint_pos_R = temp.joint_pos_R;
      
    'loaded'
    clear temp
else
    temp = pathDB3_nosave(settings,pathDB);
%     temp = load('pathDB3.mat');

    pathDB.qL1_filt1 = temp.qL1_filt1;
    pathDB.qL2_filt1 = temp.qL2_filt1;
    pathDB.qL3_filt1 = temp.qL3_filt1;
    pathDB.qL4_filt1 = temp.qL4_filt1;

    pathDB.qR1_filt1 = temp.qR1_filt1;
    pathDB.qR2_filt1 = temp.qR2_filt1;
    pathDB.qR3_filt1 = temp.qR3_filt1;
    pathDB.qR4_filt1 = temp.qR4_filt1;
    
    pathDB.b_alfa = temp.alfa_body;
    pathDB.b_beta = temp.beta_body;
    
    pathDB.b_roll = temp.roll_body;
    pathDB.b_pitch = temp.pitch_body;
    pathDB.b_yaw = temp.yaw_body;
    
    pathDB.u_body = temp.u_body;
    pathDB.v_body = temp.v_body;
    pathDB.w_body = temp.w_body;
    
    pathDB.ax_body = temp.ax_body;
    pathDB.ay_body = temp.ay_body;
    pathDB.az_body = temp.az_body;
    
    pathDB.body_l = temp.body_length;
    pathDB.wing_l = temp.wing_length;
    pathDB.joint_pos_L = temp.joint_pos_L;
    pathDB.joint_pos_R = temp.joint_pos_R;
    
%     'created + loaded'
    clear temp
end

%%

% Smooth the wing quaternions

if exist('pathDB4.mat')==2
    temp = load('pathDB4.mat');
    
    pathDB.omega1_L = temp.omega1_L;
    pathDB.omega2_L = temp.omega2_L;
    pathDB.omega3_L = temp.omega3_L;
    
    pathDB.qL1_filt2 = temp.qL1_filt2;
    pathDB.qL2_filt2 = temp.qL2_filt2;
    pathDB.qL3_filt2 = temp.qL3_filt2;
    pathDB.qL4_filt2 = temp.qL4_filt2;
    
    pathDB.omega1_R = temp.omega1_R;
    pathDB.omega2_R = temp.omega2_R;
    pathDB.omega3_R = temp.omega3_R;

    pathDB.qR1_filt2 = temp.qR1_filt2;
    pathDB.qR2_filt2 = temp.qR2_filt2;
    pathDB.qR3_filt2 = temp.qR3_filt2;
    pathDB.qR4_filt2 = temp.qR4_filt2;
      
    'loaded'
    clear temp
else
    temp = filter_wing_nosave(settings,pathDB);
%     temp = load('pathDB4.mat');
    
    pathDB.omega1_L = temp.omega1_L;
    pathDB.omega2_L = temp.omega2_L;
    pathDB.omega3_L = temp.omega3_L;

    pathDB.qL1_filt2 = temp.qL1_filt2;
    pathDB.qL2_filt2 = temp.qL2_filt2;
    pathDB.qL3_filt2 = temp.qL3_filt2;
    pathDB.qL4_filt2 = temp.qL4_filt2;
    
    pathDB.omega1_R = temp.omega1_R;
    pathDB.omega2_R = temp.omega2_R;
    pathDB.omega3_R = temp.omega3_R;

    pathDB.qR1_filt2 = temp.qR1_filt2;
    pathDB.qR2_filt2 = temp.qR2_filt2;
    pathDB.qR3_filt2 = temp.qR3_filt2;
    pathDB.qR4_filt2 = temp.qR4_filt2;
    
%     'created + loaded'
    clear temp
end

%%


% Extract wing kinematic data from the filtered quaternions.

if exist('pathDB5.mat')==2
    temp = load('pathDB5.mat');
    
    pathDB.u_wing_L = temp.u_wing_L;
    pathDB.v_wing_L = temp.v_wing_L;
    pathDB.w_wing_L = temp.w_wing_L;
    
    pathDB.u_wing_R = temp.u_wing_R;
    pathDB.v_wing_R = temp.v_wing_R;
    pathDB.w_wing_R = temp.w_wing_R;
    
    pathDB.u_wing_rot_L = temp.u_wing_rot_L;
    pathDB.v_wing_rot_L = temp.v_wing_rot_L;
    pathDB.w_wing_rot_L = temp.w_wing_rot_L;
    
    pathDB.u_wing_rot_R = temp.u_wing_rot_R;
    pathDB.v_wing_rot_R = temp.v_wing_rot_R;
    pathDB.w_wing_rot_R = temp.w_wing_rot_R;
    
    pathDB.u_wing_trans_L = temp.u_wing_trans_L;
    pathDB.v_wing_trans_L = temp.v_wing_trans_L;
    pathDB.w_wing_trans_L = temp.w_wing_trans_L;
    
    pathDB.u_wing_trans_R = temp.u_wing_trans_R;
    pathDB.v_wing_trans_R = temp.v_wing_trans_R;
    pathDB.w_wing_trans_R = temp.w_wing_trans_R;
    
    pathDB.u_wing_body_rot_L = temp.u_wing_body_rot_L;
    pathDB.v_wing_body_rot_L = temp.v_wing_body_rot_L;
    pathDB.w_wing_body_rot_L = temp.w_wing_body_rot_L;
    
    pathDB.u_wing_body_rot_R = temp.u_wing_body_rot_R;
    pathDB.v_wing_body_rot_R = temp.v_wing_body_rot_R;
    pathDB.w_wing_body_rot_R = temp.w_wing_body_rot_R;
    
    pathDB.alfa_L = temp.alfa_L;
    pathDB.beta_L = temp.beta_L;
    
    pathDB.alfa_R = temp.alfa_R;
    pathDB.beta_R = temp.beta_R;
    
    pathDB.Lwingtip = temp.Lwingtip;
    pathDB.Rwingtip = temp.Rwingtip;

    pathDB.L_wingbeat_loc = temp.Lwingbeat_loc;
    pathDB.R_wingbeat_loc = temp.Rwingbeat_loc;

    pathDB.phi_L = temp.phi_L;
    pathDB.theta_L = temp.theta_L;
    pathDB.eta_L = temp.eta_L;

    pathDB.phi_R = temp.phi_R;
    pathDB.theta_R = temp.theta_R;
    pathDB.eta_R = temp.eta_R;
      
    'loaded'
    clear temp
else
    temp = Wing_ref_frame_saveWBkin(settings,pathDB);
%     temp = load('pathDB5.mat');
    
    pathDB.u_wing_L = temp.u_wing_L;
    pathDB.v_wing_L = temp.v_wing_L;
    pathDB.w_wing_L = temp.w_wing_L;
    
    pathDB.u_wing_R = temp.u_wing_R;
    pathDB.v_wing_R = temp.v_wing_R;
    pathDB.w_wing_R = temp.w_wing_R;
    
    pathDB.u_wing_rot_L = temp.u_wing_rot_L;
    pathDB.v_wing_rot_L = temp.v_wing_rot_L;
    pathDB.w_wing_rot_L = temp.w_wing_rot_L;
    
    pathDB.u_wing_rot_R = temp.u_wing_rot_R;
    pathDB.v_wing_rot_R = temp.v_wing_rot_R;
    pathDB.w_wing_rot_R = temp.w_wing_rot_R;
    
    pathDB.u_wing_trans_L = temp.u_wing_trans_L;
    pathDB.v_wing_trans_L = temp.v_wing_trans_L;
    pathDB.w_wing_trans_L = temp.w_wing_trans_L;
    
    pathDB.u_wing_trans_R = temp.u_wing_trans_R;
    pathDB.v_wing_trans_R = temp.v_wing_trans_R;
    pathDB.w_wing_trans_R = temp.w_wing_trans_R;
    
    pathDB.u_wing_body_rot_L = temp.u_wing_body_rot_L;
    pathDB.v_wing_body_rot_L = temp.v_wing_body_rot_L;
    pathDB.w_wing_body_rot_L = temp.w_wing_body_rot_L;
    
    pathDB.u_wing_body_rot_R = temp.u_wing_body_rot_R;
    pathDB.v_wing_body_rot_R = temp.v_wing_body_rot_R;
    pathDB.w_wing_body_rot_R = temp.w_wing_body_rot_R;
    
    pathDB.alfa_L = temp.alfa_L;
    pathDB.beta_L = temp.beta_L;
    
    pathDB.alfa_R = temp.alfa_R;
    pathDB.beta_R = temp.beta_R;
    
    pathDB.Lwingtip = temp.Lwingtip;
    pathDB.Rwingtip = temp.Rwingtip;

    pathDB.L_wingbeat_loc = temp.Lwingbeat_loc;
    pathDB.R_wingbeat_loc = temp.Rwingbeat_loc;
    
    pathDB.phi_L = temp.phi_L;
    pathDB.theta_L = temp.theta_L;
    pathDB.eta_L = temp.eta_L;

    pathDB.phi_R = temp.phi_R;
    pathDB.theta_R = temp.theta_R;
    pathDB.eta_R = temp.eta_R;
    
    phi_L = temp.phi_L;
    theta_L = temp.theta_L;
    eta_L = temp.eta_L;

    phi_R = temp.phi_R;
    theta_R = temp.theta_R;
    eta_R = temp.eta_R;
    
%     'created + loaded'
    clear temp
    
end

end

%% plot wb kin
eta_L = rad2deg(eta_L) - 90;
eta_R = rad2deg(eta_R) - 90;
eta_L(eta_L<-180) = eta_L(eta_L<-180) + 360;
eta_R(eta_R<-180) = eta_R(eta_R<-180) + 360;

% figure
subplot(3,1,1)
hold off
plot(rad2deg(phi_L))
hold on
plot(rad2deg(phi_R),'r')
xlabel('frame')
ylabel('stroke angle')

subplot(3,1,2)
hold off
plot(rad2deg(theta_L))
hold on
plot(rad2deg(theta_R),'r')
xlabel('frame')
ylabel('stroke deviation')

subplot(3,1,3)
hold off
plot(eta_L)
hold on
plot(eta_R,'r')
xlabel('frame')
ylabel('wing rotation')

saveas(gca,'WBkin.fig')
saveas(gca,'WBkin.png')

