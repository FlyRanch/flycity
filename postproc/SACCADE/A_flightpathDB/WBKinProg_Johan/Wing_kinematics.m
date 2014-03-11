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

clear
close all
clc

%%

add_paths = {...
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan'; ...
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan/Kalman_Filter'; ...
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan/Plot'; ...
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan/Plot/Smoothing_plots'; ... % };
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan/Body_model'; ...
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan/Body_model/mex'; ...
             '/home/florian/Dropbox/WORK/flytracker/postproc/A_flightpathDB/WBKinProg_Johan/Type2Regression';...
             };
         
addpath(char(add_paths(1)))
addpath(char(add_paths(2)))
addpath(char(add_paths(3)))
addpath(char(add_paths(4)))
         
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/mex/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/core/');
addpath('/home/florian/Dropbox/WORK/flytracker/flytracker/results/');
addpath('/home/florian/Dropbox/WORK/flytracker/postproc');

         
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

if length(temp_dir)<3
    
    'no data files in directory'
    pause
    
else
    
    
    temp_dirnames = {temp_dir.name};
    
    % FTMmod: change folder structure
    temp_foldernames = []
    for i = 3:length(temp_dir)
         dir_now = temp_dir(i);
       
        if dir_now.isdir == 1
            cd(dir_now.name)
            seq_dir = dir('*S00*');
            
            if ~isempty(seq_dir)
                for j = 1:length(seq_dir)

                    temp_seqname = [cd,'/',seq_dir(j).name];
                    temp_foldernames = [temp_foldernames; {temp_seqname}];
                end
            end
            cd ..
        end
    end
end
    settings.sequence_names = temp_foldernames;
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
    make_pathDB1(settings,pathDB);
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
        
    'created + loaded'
    clear temp
end


%%

%Filter the body position and body orientation with a Kalman filter.

if exist('pathDB2.mat')==2
    temp = load('pathDB2.mat');
        pathDB.x_filt = temp.x;
        pathDB.y_filt = temp.y;
        pathDB.z_filt = temp.z;
        
        pathDB.u_filt = temp.u;
        pathDB.v_filt = temp.v;
        pathDB.w_filt = temp.w;

        pathDB.ax_filt = temp.ax;
        pathDB.ay_filt = temp.ay;
        pathDB.az_filt = temp.az;

        pathDB.qb1_filt = temp.q1;
        pathDB.qb2_filt = temp.q2;
        pathDB.qb3_filt = temp.q3;
        pathDB.qb4_filt = temp.q4;
        
        pathDB.b_omega1 = temp.omega_1;
        pathDB.b_omega2 = temp.omega_2;
        pathDB.b_omega3 = temp.omega_3;
        
    'loaded'
    clear temp
else
    filter_body(settings,pathDB);
    temp = load('pathDB2.mat');
        pathDB.x_filt = temp.x;
        pathDB.y_filt = temp.y;
        pathDB.z_filt = temp.z;
        
        pathDB.u_filt = temp.u;
        pathDB.v_filt = temp.v;
        pathDB.w_filt = temp.w;

        pathDB.ax_filt = temp.ax;
        pathDB.ay_filt = temp.ay;
        pathDB.az_filt = temp.az;

        pathDB.qb1_filt = temp.q1;
        pathDB.qb2_filt = temp.q2;
        pathDB.qb3_filt = temp.q3;
        pathDB.qb4_filt = temp.q4;
        
        pathDB.b_omega1 = temp.omega_1;
        pathDB.b_omega2 = temp.omega_2;
        pathDB.b_omega3 = temp.omega_3;
        
        
    'created + loaded'
    clear temp
end


%%

%Correct wing quaternions for filtered body quaternion and calculate the
%body angles psi, theta, phi, alfa and beta and their derivatives.

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
    pathDB3(settings,pathDB);
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
    
    'created + loaded'
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
    filter_wing(settings,pathDB);
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
    
    'created + loaded'
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
    Wing_ref_frame(settings,pathDB);
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
    
    'created + loaded'
    clear temp
end


%%


% Create histogram data


if exist('pathDB6.mat')==2
    
    temp = load('pathDB6.mat');
    
    % Body:
    
    pathDB.phi_body_mean = temp.phi_body_mean;
    pathDB.theta_body_mean = temp.theta_body_mean;
    pathDB.xsi_body_mean = temp.xsi_body_mean;
    
    pathDB.phi_body_sd = temp.phi_body_sd;
    pathDB.theta_body_sd = temp.theta_body_sd;
    pathDB.xsi_body_sd = temp.xsi_body_sd;
    
    pathDB.omegax_body_mean = temp.omegax_body_mean;
    pathDB.omegay_body_mean = temp.omegay_body_mean;
    pathDB.omegaz_body_mean = temp.omegaz_body_mean;
    pathDB.Omega_body_mean = temp.Omega_body_mean;
    
    pathDB.omegax_body_sd = temp.omegax_body_sd;
    pathDB.omegay_body_sd = temp.omegay_body_sd;
    pathDB.omegaz_body_sd = temp.omegaz_body_sd;
    pathDB.Omega_body_sd = temp.Omega_body_sd;
    
    pathDB.alfa_body_mean = temp.alfa_body_mean;
    pathDB.beta_body_mean = temp.beta_body_mean;
    
    pathDB.alfa_body_sd = temp.alfa_body_sd;
    pathDB.beta_body_sd = temp.beta_body_sd;
    
    pathDB.u_body_mean = temp.u_body_mean;
    pathDB.v_body_mean = temp.v_body_mean;
    pathDB.w_body_mean = temp.w_body_mean;
    pathDB.U_body_mean = temp.U_body_mean;
    
    pathDB.u_body_sd = temp.u_body_sd;
    pathDB.v_body_sd = temp.v_body_sd;
    pathDB.w_body_sd = temp.w_body_sd;
    pathDB.U_body_sd = temp.U_body_sd;
    
    pathDB.ax_body_mean = temp.ax_body_mean;
    pathDB.ay_body_mean = temp.ay_body_mean;
    pathDB.az_body_mean = temp.az_body_mean;
    pathDB.a_body_mean = temp.a_body_mean;
    
    pathDB.ax_body_sd = temp.ax_body_sd;
    pathDB.ay_body_sd = temp.ay_body_sd;
    pathDB.az_body_sd = temp.az_body_sd;
    pathDB.a_body_sd = temp.a_body_sd;
    
    % Left wing:
    
    pathDB.phi_L_up_mean = temp.phi_L_up_mean;
    pathDB.theta_L_up_mean = temp.theta_L_up_mean;
    pathDB.eta_L_up_mean = temp.eta_L_up_mean;
    
    pathDB.phi_L_up_sd = temp.phi_L_up_sd;
    pathDB.theta_L_up_sd = temp.theta_L_up_sd;
    pathDB.eta_L_up_sd = temp.eta_L_up_sd;
    
    pathDB.phi_L_down_mean = temp.phi_L_down_mean;
    pathDB.theta_L_down_mean = temp.theta_L_down_mean;
    pathDB.eta_L_down_mean = temp.eta_L_down_mean;
    
    pathDB.phi_L_down_sd = temp.phi_L_down_sd;
    pathDB.theta_L_down_sd = temp.theta_L_down_sd;
    pathDB.eta_L_down_sd = temp.eta_L_down_sd;
    
    pathDB.alfa_L_up_mean = temp.alfa_L_up_mean;
    pathDB.beta_L_up_mean = temp.beta_L_up_mean;
    
    pathDB.alfa_L_up_sd = temp.alfa_L_up_sd;
    pathDB.beta_L_up_sd = temp.beta_L_up_sd;
    
    pathDB.alfa_L_down_mean = temp.alfa_L_down_mean;
    pathDB.beta_L_down_mean = temp.beta_L_down_mean;
    
    pathDB.alfa_L_down_sd = temp.alfa_L_down_sd;
    pathDB.beta_L_down_sd = temp.beta_L_down_sd;
    
    pathDB.u_L_up_mean = temp.u_L_up_mean;
    pathDB.v_L_up_mean = temp.v_L_up_mean;
    pathDB.w_L_up_mean = temp.w_L_up_mean;
    pathDB.U_L_up_mean = temp.U_L_up_mean;
    
    pathDB.u_L_up_sd = temp.u_L_up_sd;
    pathDB.v_L_up_sd = temp.v_L_up_sd;
    pathDB.w_L_up_sd = temp.w_L_up_sd;
    pathDB.U_L_up_sd = temp.U_L_up_sd;
    
    pathDB.u_L_down_mean = temp.u_L_down_mean;
    pathDB.v_L_down_mean = temp.v_L_down_mean;
    pathDB.w_L_down_mean = temp.w_L_down_mean;
    pathDB.U_L_down_mean = temp.U_L_down_mean;
    
    pathDB.u_L_down_sd = temp.u_L_down_sd;
    pathDB.v_L_down_sd = temp.v_L_down_sd;
    pathDB.w_L_down_sd = temp.w_L_down_sd;
    pathDB.U_L_down_sd = temp.U_L_down_sd;
    
    pathDB.omegax_L_up_mean = temp.omegax_L_up_mean;
    pathDB.omegay_L_up_mean = temp.omegay_L_up_mean;
    pathDB.omegaz_L_up_mean = temp.omegaz_L_up_mean;
    pathDB.Omega_L_up_mean = temp.Omega_L_up_mean;
    
    pathDB.omegax_L_up_sd = temp.omegax_L_up_sd;
    pathDB.omegay_L_up_sd = temp.omegay_L_up_sd;
    pathDB.omegaz_L_up_sd = temp.omegaz_L_up_sd;
    pathDB.Omega_L_up_sd = temp.Omega_L_up_sd;
    
    pathDB.omegax_L_down_mean = temp.omegax_L_down_mean;
    pathDB.omegay_L_down_mean = temp.omegay_L_down_mean;
    pathDB.omegaz_L_down_mean = temp.omegaz_L_down_mean;
    pathDB.Omega_L_down_mean = temp.Omega_L_down_mean;
    
    pathDB.omegax_L_down_sd = temp.omegax_L_down_sd;
    pathDB.omegay_L_down_sd = temp.omegay_L_down_sd;
    pathDB.omegaz_L_down_sd = temp.omegaz_L_down_sd;
    pathDB.Omega_L_down_sd = temp.Omega_L_down_sd;
    
    
    % Right wing:
    
    pathDB.phi_R_up_mean = temp.phi_R_up_mean;
    pathDB.theta_R_up_mean = temp.theta_R_up_mean;
    pathDB.eta_R_up_mean = temp.eta_R_up_mean;
    
    pathDB.phi_R_up_sd = temp.phi_R_up_sd;
    pathDB.theta_R_up_sd = temp.theta_R_up_sd;
    pathDB.eta_R_up_sd = temp.eta_R_up_sd;
    
    pathDB.phi_R_down_mean = temp.phi_R_down_mean;
    pathDB.theta_R_down_mean = temp.theta_R_down_mean;
    pathDB.eta_R_down_mean = temp.eta_R_down_mean;
    
    pathDB.phi_R_down_sd = temp.phi_R_down_sd;
    pathDB.theta_R_down_sd = temp.theta_R_down_sd;
    pathDB.eta_R_down_sd = temp.eta_R_down_sd;
    
    pathDB.alfa_R_up_mean = temp.alfa_R_up_mean;
    pathDB.beta_R_up_mean = temp.beta_R_up_mean;
    
    pathDB.alfa_R_up_sd = temp.alfa_R_up_sd;
    pathDB.beta_R_up_sd = temp.beta_R_up_sd;
    
    pathDB.alfa_R_down_mean = temp.alfa_R_down_mean;
    pathDB.beta_R_down_mean = temp.beta_R_down_mean;
    
    pathDB.alfa_R_down_sd = temp.alfa_R_down_sd;
    pathDB.beta_R_down_sd = temp.beta_R_down_sd;
    
    pathDB.u_R_up_mean = temp.u_R_up_mean;
    pathDB.v_R_up_mean = temp.v_R_up_mean;
    pathDB.w_R_up_mean = temp.w_R_up_mean;
    pathDB.U_R_up_mean = temp.U_R_up_mean;
    
    pathDB.u_R_up_sd = temp.u_R_up_sd;
    pathDB.v_R_up_sd = temp.v_R_up_sd;
    pathDB.w_R_up_sd = temp.w_R_up_sd;
    pathDB.U_R_up_sd = temp.U_R_up_sd;
    
    pathDB.u_R_down_mean = temp.u_R_down_mean;
    pathDB.v_R_down_mean = temp.v_R_down_mean;
    pathDB.w_R_down_mean = temp.w_R_down_mean;
    pathDB.U_R_down_mean = temp.U_R_down_mean;
    
    pathDB.u_R_down_sd = temp.u_R_down_sd;
    pathDB.v_R_down_sd = temp.v_R_down_sd;
    pathDB.w_R_down_sd = temp.w_R_down_sd;
    pathDB.U_R_down_sd = temp.U_R_down_sd;
    
    pathDB.omegax_R_up_mean = temp.omegax_R_up_mean;
    pathDB.omegay_R_up_mean = temp.omegay_R_up_mean;
    pathDB.omegaz_R_up_mean = temp.omegaz_R_up_mean;
    pathDB.Omega_R_up_mean = temp.Omega_R_up_mean;
    
    pathDB.omegax_R_up_sd = temp.omegax_R_up_sd;
    pathDB.omegay_R_up_sd = temp.omegay_R_up_sd;
    pathDB.omegaz_R_up_sd = temp.omegaz_R_up_sd;
    pathDB.Omega_R_up_sd = temp.Omega_R_up_sd;
    
    pathDB.omegax_R_down_mean = temp.omegax_R_down_mean;
    pathDB.omegay_R_down_mean = temp.omegay_R_down_mean;
    pathDB.omegaz_R_down_mean = temp.omegaz_R_down_mean;
    pathDB.Omega_R_down_mean = temp.Omega_R_down_mean;
    
    pathDB.omegax_R_down_sd = temp.omegax_R_down_sd;
    pathDB.omegay_R_down_sd = temp.omegay_R_down_sd;
    pathDB.omegaz_R_down_sd = temp.omegaz_R_down_sd;
    pathDB.Omega_R_down_sd = temp.Omega_R_down_sd;
    
      
    'loaded'
    clear temp
else
    PathDB6(settings,pathDB);
    
    temp = load('pathDB6.mat');
    
        % Body:
    
    pathDB.phi_body_mean = temp.phi_body_mean;
    pathDB.theta_body_mean = temp.theta_body_mean;
    pathDB.xsi_body_mean = temp.xsi_body_mean;
    
    pathDB.phi_body_sd = temp.phi_body_sd;
    pathDB.theta_body_sd = temp.theta_body_sd;
    pathDB.xsi_body_sd = temp.xsi_body_sd;
    
    pathDB.omegax_body_mean = temp.omegax_body_mean;
    pathDB.omegay_body_mean = temp.omegay_body_mean;
    pathDB.omegaz_body_mean = temp.omegaz_body_mean;
    pathDB.Omega_body_mean = temp.Omega_body_mean;
    
    pathDB.omegax_body_sd = temp.omegax_body_sd;
    pathDB.omegay_body_sd = temp.omegay_body_sd;
    pathDB.omegaz_body_sd = temp.omegaz_body_sd;
    pathDB.Omega_body_sd = temp.Omega_body_sd;
    
    pathDB.alfa_body_mean = temp.alfa_body_mean;
    pathDB.beta_body_mean = temp.beta_body_mean;
    
    pathDB.alfa_body_sd = temp.alfa_body_sd;
    pathDB.beta_body_sd = temp.beta_body_sd;
    
    pathDB.u_body_mean = temp.u_body_mean;
    pathDB.v_body_mean = temp.v_body_mean;
    pathDB.w_body_mean = temp.w_body_mean;
    pathDB.U_body_mean = temp.U_body_mean;
    
    pathDB.u_body_sd = temp.u_body_sd;
    pathDB.v_body_sd = temp.v_body_sd;
    pathDB.w_body_sd = temp.w_body_sd;
    pathDB.U_body_sd = temp.U_body_sd;
    
    pathDB.ax_body_mean = temp.ax_body_mean;
    pathDB.ay_body_mean = temp.ay_body_mean;
    pathDB.az_body_mean = temp.az_body_mean;
    pathDB.a_body_mean = temp.a_body_mean;
    
    pathDB.ax_body_sd = temp.ax_body_sd;
    pathDB.ay_body_sd = temp.ay_body_sd;
    pathDB.az_body_sd = temp.az_body_sd;
    pathDB.a_body_sd = temp.a_body_sd;
    
    % Left wing:
    
    pathDB.phi_L_up_mean = temp.phi_L_up_mean;
    pathDB.theta_L_up_mean = temp.theta_L_up_mean;
    pathDB.eta_L_up_mean = temp.eta_L_up_mean;
    
    pathDB.phi_L_up_sd = temp.phi_L_up_sd;
    pathDB.theta_L_up_sd = temp.theta_L_up_sd;
    pathDB.eta_L_up_sd = temp.eta_L_up_sd;
    
    pathDB.phi_L_down_mean = temp.phi_L_down_mean;
    pathDB.theta_L_down_mean = temp.theta_L_down_mean;
    pathDB.eta_L_down_mean = temp.eta_L_down_mean;
    
    pathDB.phi_L_down_sd = temp.phi_L_down_sd;
    pathDB.theta_L_down_sd = temp.theta_L_down_sd;
    pathDB.eta_L_down_sd = temp.eta_L_down_sd;
    
    pathDB.alfa_L_up_mean = temp.alfa_L_up_mean;
    pathDB.beta_L_up_mean = temp.beta_L_up_mean;
    
    pathDB.alfa_L_up_sd = temp.alfa_L_up_sd;
    pathDB.beta_L_up_sd = temp.beta_L_up_sd;
    
    pathDB.alfa_L_down_mean = temp.alfa_L_down_mean;
    pathDB.beta_L_down_mean = temp.beta_L_down_mean;
    
    pathDB.alfa_L_down_sd = temp.alfa_L_down_sd;
    pathDB.beta_L_down_sd = temp.beta_L_down_sd;
    
    pathDB.u_L_up_mean = temp.u_L_up_mean;
    pathDB.v_L_up_mean = temp.v_L_up_mean;
    pathDB.w_L_up_mean = temp.w_L_up_mean;
    pathDB.U_L_up_mean = temp.U_L_up_mean;
    
    pathDB.u_L_up_sd = temp.u_L_up_sd;
    pathDB.v_L_up_sd = temp.v_L_up_sd;
    pathDB.w_L_up_sd = temp.w_L_up_sd;
    pathDB.U_L_up_sd = temp.U_L_up_sd;
    
    pathDB.u_L_down_mean = temp.u_L_down_mean;
    pathDB.v_L_down_mean = temp.v_L_down_mean;
    pathDB.w_L_down_mean = temp.w_L_down_mean;
    pathDB.U_L_down_mean = temp.U_L_down_mean;
    
    pathDB.u_L_down_sd = temp.u_L_down_sd;
    pathDB.v_L_down_sd = temp.v_L_down_sd;
    pathDB.w_L_down_sd = temp.w_L_down_sd;
    pathDB.U_L_down_sd = temp.U_L_down_sd;
    
    pathDB.omegax_L_up_mean = temp.omegax_L_up_mean;
    pathDB.omegay_L_up_mean = temp.omegay_L_up_mean;
    pathDB.omegaz_L_up_mean = temp.omegaz_L_up_mean;
    pathDB.Omega_L_up_mean = temp.Omega_L_up_mean;
    
    pathDB.omegax_L_up_sd = temp.omegax_L_up_sd;
    pathDB.omegay_L_up_sd = temp.omegay_L_up_sd;
    pathDB.omegaz_L_up_sd = temp.omegaz_L_up_sd;
    pathDB.Omega_L_up_sd = temp.Omega_L_up_sd;
    
    pathDB.omegax_L_down_mean = temp.omegax_L_down_mean;
    pathDB.omegay_L_down_mean = temp.omegay_L_down_mean;
    pathDB.omegaz_L_down_mean = temp.omegaz_L_down_mean;
    pathDB.Omega_L_down_mean = temp.Omega_L_down_mean;
    
    pathDB.omegax_L_down_sd = temp.omegax_L_down_sd;
    pathDB.omegay_L_down_sd = temp.omegay_L_down_sd;
    pathDB.omegaz_L_down_sd = temp.omegaz_L_down_sd;
    pathDB.Omega_L_down_sd = temp.Omega_L_down_sd;
    
    
    % Right wing:
    
    pathDB.phi_R_up_mean = temp.phi_R_up_mean;
    pathDB.theta_R_up_mean = temp.theta_R_up_mean;
    pathDB.eta_R_up_mean = temp.eta_R_up_mean;
    
    pathDB.phi_R_up_sd = temp.phi_R_up_sd;
    pathDB.theta_R_up_sd = temp.theta_R_up_sd;
    pathDB.eta_R_up_sd = temp.eta_R_up_sd;
    
    pathDB.phi_R_down_mean = temp.phi_R_down_mean;
    pathDB.theta_R_down_mean = temp.theta_R_down_mean;
    pathDB.eta_R_down_mean = temp.eta_R_down_mean;
    
    pathDB.phi_R_down_sd = temp.phi_R_down_sd;
    pathDB.theta_R_down_sd = temp.theta_R_down_sd;
    pathDB.eta_R_down_sd = temp.eta_R_down_sd;
    
    pathDB.alfa_R_up_mean = temp.alfa_R_up_mean;
    pathDB.beta_R_up_mean = temp.beta_R_up_mean;
    
    pathDB.alfa_R_up_sd = temp.alfa_R_up_sd;
    pathDB.beta_R_up_sd = temp.beta_R_up_sd;
    
    pathDB.alfa_R_down_mean = temp.alfa_R_down_mean;
    pathDB.beta_R_down_mean = temp.beta_R_down_mean;
    
    pathDB.alfa_R_down_sd = temp.alfa_R_down_sd;
    pathDB.beta_R_down_sd = temp.beta_R_down_sd;
    
    pathDB.u_R_up_mean = temp.u_R_up_mean;
    pathDB.v_R_up_mean = temp.v_R_up_mean;
    pathDB.w_R_up_mean = temp.w_R_up_mean;
    pathDB.U_R_up_mean = temp.U_R_up_mean;
    
    pathDB.u_R_up_sd = temp.u_R_up_sd;
    pathDB.v_R_up_sd = temp.v_R_up_sd;
    pathDB.w_R_up_sd = temp.w_R_up_sd;
    pathDB.U_R_up_sd = temp.U_R_up_sd;
    
    pathDB.u_R_down_mean = temp.u_R_down_mean;
    pathDB.v_R_down_mean = temp.v_R_down_mean;
    pathDB.w_R_down_mean = temp.w_R_down_mean;
    pathDB.U_R_down_mean = temp.U_R_down_mean;
    
    pathDB.u_R_down_sd = temp.u_R_down_sd;
    pathDB.v_R_down_sd = temp.v_R_down_sd;
    pathDB.w_R_down_sd = temp.w_R_down_sd;
    pathDB.U_R_down_sd = temp.U_R_down_sd;
    
    pathDB.omegax_R_up_mean = temp.omegax_R_up_mean;
    pathDB.omegay_R_up_mean = temp.omegay_R_up_mean;
    pathDB.omegaz_R_up_mean = temp.omegaz_R_up_mean;
    pathDB.Omega_R_up_mean = temp.Omega_R_up_mean;
    
    pathDB.omegax_R_up_sd = temp.omegax_R_up_sd;
    pathDB.omegay_R_up_sd = temp.omegay_R_up_sd;
    pathDB.omegaz_R_up_sd = temp.omegaz_R_up_sd;
    pathDB.Omega_R_up_sd = temp.Omega_R_up_sd;
    
    pathDB.omegax_R_down_mean = temp.omegax_R_down_mean;
    pathDB.omegay_R_down_mean = temp.omegay_R_down_mean;
    pathDB.omegaz_R_down_mean = temp.omegaz_R_down_mean;
    pathDB.Omega_R_down_mean = temp.Omega_R_down_mean;
    
    pathDB.omegax_R_down_sd = temp.omegax_R_down_sd;
    pathDB.omegay_R_down_sd = temp.omegay_R_down_sd;
    pathDB.omegaz_R_down_sd = temp.omegaz_R_down_sd;
    pathDB.Omega_R_down_sd = temp.Omega_R_down_sd;
    

    
    
    'created + loaded'
    clear temp
end





%%

% Quasi steady model:

if exist('pathDB7.mat')==2
    temp = load('pathDB7.mat');
        
    pathDB.F_joint_L_down = temp.F_joint_L_down;
    pathDB.M_joint_L_down = temp.M_joint_L_down;
    
    pathDB.F_joint_L_up = temp.F_joint_L_up;
    pathDB.M_joint_L_up = temp.M_joint_L_up;
    
    pathDB.F_joint_R_down = temp.F_joint_R_down;
    pathDB.M_joint_R_down = temp.M_joint_R_down;
    
    pathDB.F_joint_R_up = temp.F_joint_R_up;
    pathDB.M_joint_R_up = temp.M_joint_R_up;
    
    pathDB.F_cg = temp.F_cg;
    pathDB.M_cg = temp.M_cg;
    
    pathDB.R_L_down = temp.R_L_down;
    pathDB.Fn_L_down = temp.Fn_L_down;
    pathDB.Ft_L_down = temp.Ft_L_down;
    pathDB.Mn_L_down = temp.Mn_L_down;
    pathDB.Mt_L_down = temp.Mt_L_down;
    
    pathDB.R_L_up = temp.R_L_up;
    pathDB.Fn_L_up = temp.Fn_L_up;
    pathDB.Ft_L_up = temp.Ft_L_up;
    pathDB.Mn_L_up = temp.Mn_L_up;
    pathDB.Mt_L_up = temp.Mt_L_up;
    
    pathDB.R_R_down = temp.R_R_down;
    pathDB.Fn_R_down = temp.Fn_R_down;
    pathDB.Ft_R_down = temp.Ft_R_down;
    pathDB.Mn_R_down = temp.Mn_R_down;
    pathDB.Mt_R_down = temp.Mt_R_down;
    
    pathDB.R_R_up = temp.R_R_up;
    pathDB.Fn_R_up = temp.Fn_R_up;
    pathDB.Ft_R_up = temp.Ft_R_up;
    pathDB.Mn_R_up = temp.Mn_R_up;
    pathDB.Mt_R_up = temp.Mt_R_up;
    
    pathDB.F_L = temp.F_L;
    pathDB.F_R = temp.F_R;
    
    pathDB.F_L_arm = temp.F_L_arm;
    pathDB.F_R_arm = temp.F_R_arm;
      
    pathDB.m_est_point = temp.m_est_point;
    
    pathDB.m_est_dynamical_model = temp.m_est_dynamical_model;
    
    pathDB.r_body_L_down = temp.r_body_L_down;
    pathDB.r_body_L_up = temp.r_body_L_up;
    
    pathDB.r_body_R_down = temp.r_body_R_down;
    pathDB.r_body_R_up = temp.r_body_R_up;
    
    pathDB.q_avg_down_body = temp.q_avg_down_body;
    pathDB.q_avg_up_body = temp.q_avg_up_body;
    
    pathDB.q_avg_down_L = temp.q_avg_down_L;
    pathDB.q_avg_up_L = temp.q_avg_up_L;
    
    pathDB.q_avg_down_R = temp.q_avg_down_R;
    pathDB.q_avg_up_R = temp.q_avg_up_R;    
    
    pathDB.down_time_L = temp.down_time_L;
    pathDB.up_time_L = temp.up_time_L;
        
    pathDB.down_time_R = temp.down_time_R;
    pathDB.up_time_R = temp.up_time_R;
    
    pathDB.wingbeat_time = temp.wingbeat_time;
        
    'loaded'
    clear temp
else
    pathDB7(settings,pathDB);
    temp = load('pathDB7.mat');
        
    pathDB.F_joint_L_down = temp.F_joint_L_down;
    pathDB.M_joint_L_down = temp.M_joint_L_down;
    
    pathDB.F_joint_L_up = temp.F_joint_L_up;
    pathDB.M_joint_L_up = temp.M_joint_L_up;
    
    pathDB.F_joint_R_down = temp.F_joint_R_down;
    pathDB.M_joint_R_down = temp.M_joint_R_down;
    
    pathDB.F_joint_R_up = temp.F_joint_R_up;
    pathDB.M_joint_R_up = temp.M_joint_R_up;
    
    pathDB.F_cg = temp.F_cg;
    pathDB.M_cg = temp.M_cg;
    
    pathDB.R_L_down = temp.R_L_down;
    pathDB.Fn_L_down = temp.Fn_L_down;
    pathDB.Ft_L_down = temp.Ft_L_down;
    pathDB.Mn_L_down = temp.Mn_L_down;
    pathDB.Mt_L_down = temp.Mt_L_down;
    
    pathDB.R_L_up = temp.R_L_up;
    pathDB.Fn_L_up = temp.Fn_L_up;
    pathDB.Ft_L_up = temp.Ft_L_up;
    pathDB.Mn_L_up = temp.Mn_L_up;
    pathDB.Mt_L_up = temp.Mt_L_up;
    
    pathDB.R_R_down = temp.R_R_down;
    pathDB.Fn_R_down = temp.Fn_R_down;
    pathDB.Ft_R_down = temp.Ft_R_down;
    pathDB.Mn_R_down = temp.Mn_R_down;
    pathDB.Mt_R_down = temp.Mt_R_down;
    
    pathDB.R_R_up = temp.R_R_up;
    pathDB.Fn_R_up = temp.Fn_R_up;
    pathDB.Ft_R_up = temp.Ft_R_up;
    pathDB.Mn_R_up = temp.Mn_R_up;
    pathDB.Mt_R_up = temp.Mt_R_up;
    
    pathDB.F_L = temp.F_L;
    pathDB.F_R = temp.F_R;
    
    pathDB.F_L_arm = temp.F_L_arm;
    pathDB.F_R_arm = temp.F_R_arm;
    
    pathDB.m_est_point = temp.m_est_point;
    
    pathDB.m_est_dynamical_model = temp.m_est_dynamical_model;
    
    pathDB.r_body_L_down = temp.r_body_L_down;
    pathDB.r_body_L_up = temp.r_body_L_up;
    
    pathDB.r_body_R_down = temp.r_body_R_down;
    pathDB.r_body_R_up = temp.r_body_R_up;
    
    pathDB.q_avg_down_body = temp.q_avg_down_body;
    pathDB.q_avg_up_body = temp.q_avg_up_body;
    
    pathDB.q_avg_down_L = temp.q_avg_down_L;
    pathDB.q_avg_up_L = temp.q_avg_up_L;
    
    pathDB.q_avg_down_R = temp.q_avg_down_R;
    pathDB.q_avg_up_R = temp.q_avg_up_R;
    
    pathDB.down_time_L = temp.down_time_L;
    pathDB.up_time_L = temp.up_time_L;
        
    pathDB.down_time_R = temp.down_time_R;
    pathDB.up_time_R = temp.up_time_R;
    
    pathDB.wingbeat_time = temp.wingbeat_time;
        
    'created + loaded';
    clear temp
end


%%

wing_kinematics_savedata(settings,pathDB)

%%
% Dynamic_model_fruit_fly4(settings, pathDB, 29)

%Dynamic_model_fruit_fly3(settings, pathDB, 13)

% Fitting_Florian( settings, pathDB)

% Wingkinematics_Legendre(settings, pathDB)

% Dynamic_test(settings, pathDB)

% % for j = 1:size(pathDB.x,2)
% 
% seq_nr = 1
% 
% 
% [F_r, M_r, T_b_trans, T_b_rot, T_L_trans, T_L_rot, T_R_trans, T_R_rot, F_r_avg, M_r_avg, t_avg, M_r_meas, M_r_meas_avg, F_r_global] = Dynamic_model_fruit_fly(settings, pathDB, seq_nr);
% 
% 
% 
% % figure()
% % plot(F_r(1,:),'Color','b')
% % hold on
% % plot(F_r(2,:),'Color','r')
% % plot(F_r(3,:),'Color','g')
% % legend('Fx','Fy','Fz')
% % hold off
% 
% 
% wing_l = pathDB.wing_l(seq_nr);
%     
% m_wl = 1e-6*0.1078*wing_l^3.008;
% 
% start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
% stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
% 
% t = pathDB.t(start:stop);
% 
% dt = t(2)-t(1);
% 
% ax = pathDB.ax_body(start:stop,seq_nr);
% ay = pathDB.ay_body(start:stop,seq_nr);
% az = pathDB.az_body(start:stop,seq_nr);
% 
% % figure()
% % subplot(3,1,1); plot(t,1e3*F_r(1,:)./m_wl,t,ax,pathDB.t(start+t_avg),1e3*F_r_avg(1,:)./m_wl,'o')
% % subplot(3,1,2); plot(t,1e3*F_r(2,:)./m_wl,t,ay,pathDB.t(start+t_avg),1e3*F_r_avg(2,:)./m_wl,'o')
% % subplot(3,1,3); plot(t,1e3*F_r(3,:)./m_wl,t,az,pathDB.t(start+t_avg),1e3*F_r_avg(3,:)./m_wl,'o')
% 
% % figure()
% % plot(M_r(1,:),'Color','b')
% % hold on
% % plot(M_r(2,:),'Color','r')
% % plot(M_r(3,:),'Color','g')
% % legend('Mx','My','Mz')
% % hold off
% 
% figure()
% subplot(3,1,1); plot(t,0,t,1e-3.*ax.*m_wl,pathDB.t(start+t_avg),F_r_avg(1,:))
% subplot(3,1,2); plot(t,0,t,1e-3.*ay.*m_wl,pathDB.t(start+t_avg),F_r_avg(2,:))
% subplot(3,1,3); plot(t,0,t,1e-3.*az.*m_wl,pathDB.t(start+t_avg),F_r_avg(3,:))
% 
% figure()
% subplot(3,1,1); plot(t,0,t,ax,pathDB.t(start+t_avg),1e3*F_r_avg(1,:)./m_wl)
% subplot(3,1,2); plot(t,0,t,ay,pathDB.t(start+t_avg),1e3*F_r_avg(2,:)./m_wl)
% subplot(3,1,3); plot(t,0,t,az,pathDB.t(start+t_avg),1e3*F_r_avg(3,:)./m_wl)
% 
% 
% figure()
% subplot(3,1,1); plot(t,M_r(1,:),pathDB.t(start+t_avg),M_r_avg(1,:),t,M_r_meas(1,:),pathDB.t(start+t_avg),M_r_meas_avg(1,:))
% subplot(3,1,2); plot(t,M_r(2,:),pathDB.t(start+t_avg),M_r_avg(2,:),t,M_r_meas(2,:),pathDB.t(start+t_avg),M_r_meas_avg(2,:))
% subplot(3,1,3); plot(t,M_r(3,:),pathDB.t(start+t_avg),M_r_avg(3,:),t,M_r_meas(3,:),pathDB.t(start+t_avg),M_r_meas_avg(3,:))
% 
% figure()
% subplot(3,1,1); plot(pathDB.t(start+t_avg),M_r_avg(1,:),pathDB.t(start+t_avg),M_r_meas_avg(1,:))
% subplot(3,1,2); plot(pathDB.t(start+t_avg),M_r_avg(2,:),pathDB.t(start+t_avg),M_r_meas_avg(2,:))
% subplot(3,1,3); plot(pathDB.t(start+t_avg),M_r_avg(3,:),pathDB.t(start+t_avg),M_r_meas_avg(3,:))
% 
% figure()
% plot(T_b_trans,'Color','r')
% hold on
% plot(T_b_rot,'Color','b')
% plot(T_b_trans+T_b_rot,'Color','g')
% title('Kinetic energies body')
% legend('trans','rot','combined')
% hold off
% 
% figure()
% plot(T_L_trans,'Color','r')
% hold on
% plot(T_L_rot,'Color','b')
% plot(T_L_trans+T_L_rot,'Color','g')
% title('Kinetic energies left wing')
% legend('trans','rot','combined')
% hold off
% 
% figure()
% plot(T_R_trans,'Color','r')
% hold on
% plot(T_R_rot,'Color','b')
% plot(T_R_trans+T_R_rot,'Color','g')
% title('Kinetic energies right wing')
% legend('trans','rot','combined')
% hold off
% 
% figure()
% plot(T_b_trans+T_L_trans+T_R_trans,'Color','r')
% hold on
% plot(T_b_rot+T_L_rot+T_R_rot,'Color','b')
% plot(T_b_trans+T_L_trans+T_R_trans+T_b_rot+T_L_rot+T_R_rot,'Color','g')
% title('Kinetic energy fruit fly')
% legend('trans','rot','combined')
% hold off
% 
% figure()
% plot(T_L_rot,'Color','r')
% hold on
% plot(T_R_rot,'Color','b')
% hold off
% 
% figure()
% plot(T_L_trans,'Color','r')
% hold on
% plot(T_R_trans,'Color','b')
% hold off
% 
% 
% figure()
% plot(T_L_trans+T_R_trans)
% 
% figure()
% plot3(pathDB.x_filt(start:stop,seq_nr),pathDB.y_filt(start:stop,seq_nr),pathDB.z_filt(start:stop,seq_nr))
% hold on
% quiver3(pathDB.x_filt(start+t_avg,seq_nr),pathDB.y_filt(start+t_avg,seq_nr),pathDB.z_filt(start+t_avg,seq_nr),F_r_global(1,:)',F_r_global(2,:)',F_r_global(3,:)')
% hold off

% end


%PLOT_generator2( settings, pathDB )

% 
% ProjectFlorian(settings, pathDB)

% Return_xh(settings, pathDB)
% 
% 
% Histogram(settings, pathDB)
% 
% 
% Plot_relations(settings, pathDB)


% PLOT_generator( settings, pathDB )


% for i = 1:1:size(pathDB.x,2)
% 
% Fly_movie(settings,pathDB,i)
% 
% end

% 
% for i = 1:1:size(pathDB.x,2)
% 
% Fly_movie4(settings,pathDB,i)
% 
% end
% 
% for i = 1:1:size(pathDB.x,2)
% 
% Fly_movie5(settings,pathDB,i)
% 
% end

% Histogram_Floris(settings, pathDB)

