%% Run variables
%date_str = '20140402'
for i= 2:12
    date_str = '20140429'
    seq_number = i;
    %startframe = 193;
    BGframe = 1;
    endframe = 1365;

    %% Deployment variables
    software_root_dir = '/home/matt/flycity/flytracker_tethered'
    data_root_dir = '/media/thad01b/'

    %% Kine settings
    kine_dir = fullfile(software_root_dir,'kine');
    setup_path = fullfile(software_root_dir,'kine','setup');
    kine_mfile_path = kine_dir;
    setup_file = 'flightmuscles.mat';


    %% Flytracker settings
    % PAR is a static structure that defines various parameters of the tracking
    % sequence
    % String formating conventions
    seq_format_str = 'C%03dH%03dS%04d' %cam number, head number, sequence number
    sol_format_str = '%s_S%04d'%date_str, sol number
    % File path formating
    MoveModelPath = fullfile(software_root_dir,'Models',filesep);
    sol_root_path = fullfile(data_root_dir,'solutions');
    seq_sol_path = fullfile(sol_root_path,sprintf(sol_format_str,date_str,seq_number),filesep)

    % software paths
    addpath(fullfile(software_root_dir,'flytracker_library'));
    addpath(fullfile(software_root_dir,'mex'));
    addpath(fullfile(software_root_dir,'core'));
    addpath(fullfile(software_root_dir,'results'));
    addpath(fullfile(software_root_dir,'SpinConv'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Kalman_Filter'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Plot'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Plot','Smoothing_plots'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Body_model'));
    %addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Body_model','mex'));
    addpath(fullfile(software_root_dir,'postproc','AA_check_wbkin','WBKinProg_Johan','Type2Regression'));

    cd(seq_sol_path);
    %paste_top_projection_makemovie
    A_Wing_kinematics_v2_checkFlyData
end
