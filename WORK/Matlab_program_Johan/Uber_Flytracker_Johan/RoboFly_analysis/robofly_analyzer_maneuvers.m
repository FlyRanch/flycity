function robofly_analyzer_maneuvers( settings, pathDB )


close all;

%   Program to analyze maneuvering Robofly data.

% folder_names = {    'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_FX_back', ...
%                     'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_FX_forward', ...
%                     'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_FY', ...
%                     'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_FZ_up', ...
%                     'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_MX', ...
%                     'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_MY_up', ...
%                     'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06252013_F_MZ'};
%                 
% figure_names = { 'FXback', ...
%                  'FXforward', ...
%                  'FY', ...
%                  'FZup', ...
%                  'MX', ...
%                  'MYup', ...
%                  'MZ' };


folder_names = {    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_FX_back', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_FX_forward', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_FY', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_FZ_up', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_MX', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_MY_down', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_MY_up', ...
                    'C:/Users/Johan/Documents/Thesis_data/Johan_RoboFly_tests_29_10/F_MZ'};
                
figure_names = { 'FXback', ...
                 'FXforward', ...
                 'FY', ...
                 'FZup', ...
                 'MX', ...
                 'MYdown', ...
                 'MYup', ...
                 'MZ' };
             
plot_loc = 'C:/Users/Johan/Documents/Thesis_data/Plots';
                
                    %'C:/Users/Johan/Documents/Thesis_data/robofly_data_8_7_2013/Johan_06252013/Johan_06272013_FlowVis'    };
                    
[ ~, ~, ~, ~, ~, ~ ] = g_calib( 0, 0, 0, 0, 0, 0, 0 );


[ FM_rot ] = quasi_steady_man_wingkin( settings, pathDB, 1 );
    
FM_qr_ax_forward     = FM_rot.ax_forward;
FM_qr_ax_back        = FM_rot.ax_back;
FM_qr_ay             = FM_rot.ay;
FM_qr_az_up          = FM_rot.az_up;
FM_qr_wx             = FM_rot.wx;
FM_qr_wy_down        = FM_rot.wy_down;
FM_qr_wy_up          = FM_rot.wy_up;
FM_qr_wz             = FM_rot.wz;

FM_qr_glob           = FM_rot.glob;

[ FM_trans ] = quasi_steady_man_wingkin( settings, pathDB, 0 );
    
FM_qt_ax_forward     = FM_trans.ax_forward;
FM_qt_ax_back        = FM_trans.ax_back;
FM_qt_ay             = FM_trans.ay;
FM_qt_az_up          = FM_trans.az_up;
FM_qt_wx             = FM_trans.wx;
FM_qt_wy_down        = FM_trans.wy_down;
FM_qt_wy_up          = FM_trans.wy_up;
FM_qt_wz             = FM_trans.wz;

FM_qt_glob           = FM_trans.glob;
                
% Analyze the different folders with experiments:

for i = 1:length(folder_names)
    
   
    % Open folder i:
    
    cd(char(folder_names(i)))
    
    temp_dir = dir;
    
    dir_names = {temp_dir.name};

    temp_isdir = [temp_dir.isdir];
    
    exp_names1 = [];

    for j = 1:length(dir_names)

        if temp_isdir(j) == 0

           exp_names1 = [exp_names1; dir_names(j)];

        end

    end
    
    exp_names = {};
    
    temp_digits = strfind(char(exp_names1(1)),'_');
    
    start_digit = temp_digits(end)+1;
        
    for k = 1:length(exp_names1)
        
        temp_name = char(exp_names1(k));
        
        file_index = str2num(temp_name(start_digit:(end-4)))+1;
        
        exp_names(file_index) = exp_names1(k);
        
    end
    
    folder_names(i)
    
    Test = [];

    for m = 1:length(exp_names)

        temp_file = load((char(exp_names(m))));

        Test.([ 'test_' int2str(m) ]) = temp_file;

    end

    Test_names = fieldnames(Test);
    


    
    time_s            = nan(80,12000);
    FX_s              = nan(80,12000);
    FY_s              = nan(80,12000);
    FZ_s              = nan(80,12000);
    MX_s              = nan(80,12000);
    MY_s              = nan(80,12000);
    MZ_s              = nan(80,12000);
    theta_L_s         = nan(80,12000);
    eta_L_s           = nan(80,12000);
    phi_L_s           = nan(80,12000);
    theta_R_s         = nan(80,12000);
    eta_R_s           = nan(80,12000);
    phi_R_s           = nan(80,12000);
    freq              = nan(80,7);
    wing_length       = nan(80,7);    
    
    FX_g            = nan(80,12000);
    FY_g            = nan(80,12000);
    FZ_g            = nan(80,12000);
    MX_g            = nan(80,12000);
    MY_g            = nan(80,12000);
    MZ_g            = nan(80,12000);
    
    
    
    
    for n = 1:80
        
        tic
        
        time            = nan(1,120000);
        FX              = nan(1,120000);
        FY              = nan(1,120000);
        FZ              = nan(1,120000);
        MX              = nan(1,120000);
        MY              = nan(1,120000);
        MZ              = nan(1,120000);
        theta_L         = nan(1,120000);
        eta_L           = nan(1,120000);
        phi_L           = nan(1,120000);
        theta_R         = nan(1,120000);
        eta_R           = nan(1,120000);
        phi_R           = nan(1,120000);
        
        test_end = length(Test.(char(Test_names(n))).t);
        
        time(1:test_end)            = Test.(char(Test_names(n))).t;
        FX(1:test_end)              = Test.(char(Test_names(n))).ft(:,1);
        FY(1:test_end)              = Test.(char(Test_names(n))).ft(:,3);
        FZ(1:test_end)              = Test.(char(Test_names(n))).ft(:,2);
        MX(1:test_end)              = Test.(char(Test_names(n))).ft(:,4);
        MY(1:test_end)              = Test.(char(Test_names(n))).ft(:,6);
        MZ(1:test_end)              = Test.(char(Test_names(n))).ft(:,5);
        theta_L(1:test_end)         = Test.(char(Test_names(n))).kine(:,3);
        eta_L(1:test_end)           = Test.(char(Test_names(n))).kine(:,1);
        phi_L(1:test_end)           = Test.(char(Test_names(n))).kine(:,5);
        theta_R(1:test_end)         = Test.(char(Test_names(n))).kine(:,2);
        eta_R(1:test_end)           = Test.(char(Test_names(n))).kine(:,4);
        phi_R(1:test_end)           = Test.(char(Test_names(n))).kine(:,6);
        freq(n,:)                   = Test.(char(Test_names(n))).freq;
        wing_length(n,:)            = Test.(char(Test_names(n))).wing_length;
        
        filt_l = 50;
        
        time_filt                   = time(1:test_end);
        FX_filt                     = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,FX(1:test_end))));
        FY_filt                     = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,FY(1:test_end))));
        FZ_filt                     = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,FZ(1:test_end))));
        MX_filt                     = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,MX(1:test_end))));
        MY_filt                     = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,MY(1:test_end))));
        MZ_filt                     = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,MZ(1:test_end))));
        theta_L_filt                = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,theta_L(1:test_end))));
        eta_L_filt                  = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,eta_L(1:test_end))));
        phi_L_filt                  = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,phi_L(1:test_end))));
        theta_R_filt                = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,theta_R(1:test_end))));
        eta_R_filt                  = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,eta_R(1:test_end))));
        phi_R_filt                  = filter(ones(1,filt_l)/filt_l,1,flipud(filter(ones(1,filt_l)/filt_l,1,phi_R(1:test_end))));

        test_end_10 = length(filt_l:10:test_end);

        time_s(n,1:test_end_10)            = time(1:10:(test_end-filt_l))./((1/7)*time(test_end));
        FX_s(n,1:test_end_10)              = FX_filt(filt_l:10:test_end);
        FY_s(n,1:test_end_10)              = FY_filt(filt_l:10:test_end);
        FZ_s(n,1:test_end_10)              = FZ_filt(filt_l:10:test_end);
        MX_s(n,1:test_end_10)              = MX_filt(filt_l:10:test_end);
        MY_s(n,1:test_end_10)              = MY_filt(filt_l:10:test_end);
        MZ_s(n,1:test_end_10)              = MZ_filt(filt_l:10:test_end);
        theta_L_s(n,1:test_end_10)         = theta_L_filt(filt_l:10:test_end);
        eta_L_s(n,1:test_end_10)           = eta_L_filt(filt_l:10:test_end);
        phi_L_s(n,1:test_end_10)           = phi_L_filt(filt_l:10:test_end);
        theta_R_s(n,1:test_end_10)         = theta_R_filt(filt_l:10:test_end);
        eta_R_s(n,1:test_end_10)           = eta_R_filt(filt_l:10:test_end);
        phi_R_s(n,1:test_end_10)           = phi_R_filt(filt_l:10:test_end);

        [FX_g(n,1:test_end_10),FY_g(n,1:test_end_10),FZ_g(n,1:test_end_10),MX_g(n,1:test_end_10),MY_g(n,1:test_end_10),MZ_g(n,1:test_end_10)] = g_calib(theta_L_s(n,1:test_end_10),eta_L_s(n,1:test_end_10),phi_L_s(n,1:test_end_10),theta_R_s(n,1:test_end_10),eta_R_s(n,1:test_end_10),phi_R_s(n,1:test_end_10),1);
        
        clear time FX FY FZ MX MY MZ theta_L eta_L phi_L theta_R eta_R phi_R 
        clear time_filt FX_filt FY_filt FZ_filt MX_filt MY_filt MZ_filt theta_L_filt eta_L_filt phi_L_filt theta_R_filt eta_R_filt phi_R_filt
        
    end
    
    t_glob      = 0:(7/13993):7;
    
    Fx_glob_rot = [FM_qr_glob(1,:) FM_qr_glob(1,2:2000) FM_qr_glob(1,2:2000) FM_qr_glob(1,2:2000) FM_qr_glob(1,2:2000) FM_qr_glob(1,2:2000) FM_qr_glob(1,2:2000)];
    Fy_glob_rot = [FM_qr_glob(2,:) FM_qr_glob(2,2:2000) FM_qr_glob(2,2:2000) FM_qr_glob(2,2:2000) FM_qr_glob(2,2:2000) FM_qr_glob(2,2:2000) FM_qr_glob(2,2:2000)];
    Fz_glob_rot = [FM_qr_glob(3,:) FM_qr_glob(3,2:2000) FM_qr_glob(3,2:2000) FM_qr_glob(3,2:2000) FM_qr_glob(3,2:2000) FM_qr_glob(3,2:2000) FM_qr_glob(3,2:2000)];
    Mx_glob_rot = [FM_qr_glob(4,:) FM_qr_glob(4,2:2000) FM_qr_glob(4,2:2000) FM_qr_glob(4,2:2000) FM_qr_glob(4,2:2000) FM_qr_glob(4,2:2000) FM_qr_glob(4,2:2000)];
    My_glob_rot = [FM_qr_glob(5,:) FM_qr_glob(5,2:2000) FM_qr_glob(5,2:2000) FM_qr_glob(5,2:2000) FM_qr_glob(5,2:2000) FM_qr_glob(5,2:2000) FM_qr_glob(5,2:2000)];
    Mz_glob_rot = [FM_qr_glob(6,:) FM_qr_glob(6,2:2000) FM_qr_glob(6,2:2000) FM_qr_glob(6,2:2000) FM_qr_glob(6,2:2000) FM_qr_glob(6,2:2000) FM_qr_glob(6,2:2000)];
    
    Fx_glob_trans = [FM_qt_glob(1,:) FM_qt_glob(1,2:2000) FM_qt_glob(1,2:2000) FM_qt_glob(1,2:2000) FM_qt_glob(1,2:2000) FM_qt_glob(1,2:2000) FM_qt_glob(1,2:2000)];
    Fy_glob_trans = [FM_qt_glob(2,:) FM_qt_glob(2,2:2000) FM_qt_glob(2,2:2000) FM_qt_glob(2,2:2000) FM_qt_glob(2,2:2000) FM_qt_glob(2,2:2000) FM_qt_glob(2,2:2000)];
    Fz_glob_trans = [FM_qt_glob(3,:) FM_qt_glob(3,2:2000) FM_qt_glob(3,2:2000) FM_qt_glob(3,2:2000) FM_qt_glob(3,2:2000) FM_qt_glob(3,2:2000) FM_qt_glob(3,2:2000)];
    Mx_glob_trans = [FM_qt_glob(4,:) FM_qt_glob(4,2:2000) FM_qt_glob(4,2:2000) FM_qt_glob(4,2:2000) FM_qt_glob(4,2:2000) FM_qt_glob(4,2:2000) FM_qt_glob(4,2:2000)];
    My_glob_trans = [FM_qt_glob(5,:) FM_qt_glob(5,2:2000) FM_qt_glob(5,2:2000) FM_qt_glob(5,2:2000) FM_qt_glob(5,2:2000) FM_qt_glob(5,2:2000) FM_qt_glob(5,2:2000)];
    Mz_glob_trans = [FM_qt_glob(6,:) FM_qt_glob(6,2:2000) FM_qt_glob(6,2:2000) FM_qt_glob(6,2:2000) FM_qt_glob(6,2:2000) FM_qt_glob(6,2:2000) FM_qt_glob(6,2:2000)];
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    subplot(3,2,1); hold on
    plot(time_s(1,:),(FX_s(1,:)-FX_g(1,:))*2.2429e-5,'Color',[0.5 0.5 0.5])
    plot(t_glob,Fx_glob_rot,'r')
    plot(t_glob,Fx_glob_trans,'b')
    hold off
    xlabel('t/T')
    ylabel('Fx [N]')
    xlim([5 7])
    subplot(3,2,3); hold on
    plot(time_s(1,:),(FY_s(1,:)-FY_g(1,:))*2.2429e-5,'Color',[0.5 0.5 0.5])
    plot(t_glob,Fy_glob_rot,'r')
    plot(t_glob,Fy_glob_trans,'b')
    hold off
    xlabel('t/T')
    ylabel('Fy [N]')
    xlim([5 7])
    subplot(3,2,5); hold on
    plot(time_s(1,:),(-(FZ_s(1,:)-FZ_g(1,:)))*2.2429e-5,'Color',[0.5 0.5 0.5])
    plot(t_glob,Fz_glob_rot,'r')
    plot(t_glob,Fz_glob_trans,'b')
    hold off
    xlabel('t/T')
    ylabel('Fz [N]')
    xlim([5 7])
    subplot(3,2,2); hold on
    plot(time_s(1,:),(MX_s(1,:)-MX_g(1,:))*2.2429e-5*(3/230)*1000,'Color',[0.5 0.5 0.5])
    plot(t_glob,Mx_glob_rot,'r')
    plot(t_glob,Mx_glob_trans,'b')
    hold off
    xlabel('t/T')
    ylabel('Mx [N*mm]')
    xlim([5 7])
    subplot(3,2,4); hold on
    plot(time_s(1,:),(MY_s(1,:)-MY_g(1,:))*2.2429e-5*(3/230)*1000,'Color',[0.5 0.5 0.5])
    plot(t_glob,My_glob_rot,'r')
    plot(t_glob,My_glob_trans,'b')
    hold off
    xlabel('t/T')
    ylabel('My [N*mm]')
    xlim([5 7])
    subplot(3,2,6); hold on
    plot(time_s(1,:),(MZ_s(1,:)-MZ_g(1,:))*2.2429e-5*(3/230)*1000,'Color',[0.5 0.5 0.5])
    plot(t_glob,Mz_glob_rot,'r')
    plot(t_glob,Mz_glob_trans,'b')
    hold off
    xlabel('t/T')
    ylabel('Mz [N*mm]')
    xlim([5 7])
    hold off
    legend('RoboFly','rot','trans')
    [~,h1] = suplabel('Forces and Moments average wingbeat', 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'FM_glob.eps'] ,'-depsc2');
    saveas(1, [char(plot_loc) '/RoboFly_results/' 'FM_glob' ], 'fig')    
    
    pause
        
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 1:20
        hold on
        subplot(3,2,1); plot(time_s(w,:),FX_s(w,:)-FX_g(w,:),'Color',[ (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fx [N]')
    xlim([5 7])
    ylim([-2 3])
    for w = 1:20
        hold on
        subplot(3,2,3); plot(time_s(w,:),FY_s(w,:)-FY_g(w,:),'Color',[ (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fy [N]')
    xlim([5 7])
    ylim([-1.5 1.5])
    for w = 1:20
        hold on
        subplot(3,2,5); plot(time_s(w,:),-(FZ_s(w,:)-FZ_g(w,:)),'Color',[ (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fz [N]')
    xlim([5 7])
    ylim([-3.5 1])
    for w = 1:20
        hold on
        subplot(3,2,2); plot(time_s(w,:),MX_s(w,:)-MX_g(w,:),'Color',[ (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mx [N*m]')
    xlim([5 7])
    ylim([-0.15 0.25])
    for w = 1:20
        hold on
        subplot(3,2,4); plot(time_s(w,:),MY_s(w,:)-MY_g(w,:),'Color',[ (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('My [N*m]')
    xlim([5 7])
    ylim([-0.25 0.1])
    for w = 1:20
        hold on
        subplot(3,2,6); plot(time_s(w,:),MZ_s(w,:)-MZ_g(w,:),'Color',[ (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mz [N*m]')
    xlim([5 7])
    ylim([-0.25 0.2])
    hold off
    [~,h1] = suplabel(['Forces and Moments ' char(figure_names(i))], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'FM_complete_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(1, [char(plot_loc) '/RoboFly_results/' 'FM_complete_' char(figure_names(i))], 'fig')

    figure(2)
    hFig = figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 1:20
        hold on
        subplot(3,2,1); plot(time_s(w,:),radtodeg(theta_L_s(w,:)),'Color',[ (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_L [deg]')
    for w = 1:20
        hold on
        subplot(3,2,3); plot(time_s(w,:),90-radtodeg(eta_L_s(w,:)),'Color',[ (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_L [deg]')
    for w = 1:20
        hold on
        subplot(3,2,5); plot(time_s(w,:),-radtodeg(phi_L_s(w,:)),'Color',[ (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_L [deg]')
    for w = 1:20
        hold on
        subplot(3,2,2); plot(time_s(w,:),radtodeg(theta_R_s(w,:)),'Color',[ (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_R [deg]')
    for w = 1:20
        hold on
        subplot(3,2,4); plot(time_s(w,:),90-radtodeg(eta_R_s(w,:)),'Color',[ (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_R [deg]')
    for w = 1:20
        hold on
        subplot(3,2,6); plot(time_s(w,:),-radtodeg(phi_R_s(w,:)),'Color',[ (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_R [deg]')
    hold off
    [~,h1] = suplabel(['Wing kinematics ' char(figure_names(i))], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'Wingkin_complete_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(2, [char(plot_loc) '/RoboFly_results/' 'Wingkin_complete_' char(figure_names(i))], 'fig')
    
    
    figure(3)
    hFig = figure(3);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 21:40
        hold on
        subplot(3,2,1); plot(time_s(w,:),FX_s(w,:)-FX_g(w,:),'Color',[ (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fx [N]')
    xlim([5 7])
    ylim([-2 3])
    for w = 21:40
        hold on
        subplot(3,2,3); plot(time_s(w,:),FY_s(w,:)-FY_g(w,:),'Color',[ (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fy [N]')
    xlim([5 7])
    ylim([-1.5 1.5])
    for w = 21:40
        hold on
        subplot(3,2,5); plot(time_s(w,:),-(FZ_s(w,:)-FZ_g(w,:)),'Color',[ (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fz [N]')
    xlim([5 7])
    ylim([-3.5 1])
    for w = 21:40
        hold on
        subplot(3,2,2); plot(time_s(w,:),MX_s(w,:)-MX_g(w,:),'Color',[ (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mx [N*m]')
    xlim([5 7])
    ylim([-0.15 0.25])
    for w = 21:40
        hold on
        subplot(3,2,4); plot(time_s(w,:),MY_s(w,:)-MY_g(w,:),'Color',[ (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('My [N*m]')
    xlim([5 7])
    ylim([-0.25 0.1])
    for w = 21:40
        hold on
        subplot(3,2,6); plot(time_s(w,:),MZ_s(w,:)-MZ_g(w,:),'Color',[ (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mz [N*m]')
    xlim([5 7])
    ylim([-0.25 0.2])
    hold off
    [~,h1] = suplabel(['Forces and Moments ' char(figure_names(i)) ', \theta only'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'FM_theta_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(3, [char(plot_loc) '/RoboFly_results/' 'FM_theta_' char(figure_names(i))], 'fig')

    figure(4)
    hFig = figure(4);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 21:40
        hold on
        subplot(3,2,1); plot(time_s(w,:),radtodeg(theta_L_s(w,:)),'Color',[ (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_L [deg]')
    for w = 21:40
        hold on
        subplot(3,2,3); plot(time_s(w,:),90-radtodeg(eta_L_s(w,:)),'Color',[ (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_L [deg]')
    for w = 21:40
        hold on
        subplot(3,2,5); plot(time_s(w,:),-radtodeg(phi_L_s(w,:)),'Color',[ (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_L [deg]')
    for w = 21:40
        hold on
        subplot(3,2,2); plot(time_s(w,:),radtodeg(theta_R_s(w,:)),'Color',[ (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_R [deg]')
    for w = 21:40
        hold on
        subplot(3,2,4); plot(time_s(w,:),90-radtodeg(eta_R_s(w,:)),'Color',[ (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_R [deg]')
    for w = 21:40
        hold on
        subplot(3,2,6); plot(time_s(w,:),-radtodeg(phi_R_s(w,:)),'Color',[ (0.5-(0.5*(w-20))/21) (0.5-(0.5*(w-20))/21) (0.5+(0.5*(w-20))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_R [deg]')
    hold off
    [~,h1] = suplabel(['Wing kinematics ' char(figure_names(i)) ', \theta only'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'Wingkin_theta_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(4, [char(plot_loc) '/RoboFly_results/' 'Wingkin_theta_' char(figure_names(i))], 'fig')
    

    figure(5)
    hFig = figure(5);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 41:60
        hold on
        subplot(3,2,1); plot(time_s(w,:),FX_s(w,:)-FX_g(w,:),'Color',[ (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fx [N]')
    xlim([5 7])
    ylim([-2 3])
    for w = 41:60
        hold on
        subplot(3,2,3); plot(time_s(w,:),FY_s(w,:)-FY_g(w,:),'Color',[ (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fy [N]')
    xlim([5 7])
    ylim([-1.5 1.5])
    for w = 41:60
        hold on
        subplot(3,2,5); plot(time_s(w,:),-(FZ_s(w,:)-FZ_g(w,:)),'Color',[ (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fz [N]')
    xlim([5 7])
    ylim([-3.5 1])
    for w = 41:60
        hold on
        subplot(3,2,2); plot(time_s(w,:),MX_s(w,:)-MX_g(w,:),'Color',[ (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mx [N*m]')
    xlim([5 7])
    ylim([-0.15 0.25])
    for w = 41:60
        hold on
        subplot(3,2,4); plot(time_s(w,:),MY_s(w,:)-MY_g(w,:),'Color',[ (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('My [N*m]')
    xlim([5 7])
    ylim([-0.25 0.1])
    for w = 41:60
        hold on
        subplot(3,2,6); plot(time_s(w,:),MZ_s(w,:)-MZ_g(w,:),'Color',[ (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mz [N*m]')
    xlim([5 7])
    ylim([-0.25 0.2])
    hold off
    [~,h1] = suplabel(['Forces and Moments ' char(figure_names(i)) ', \eta only'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'FM_eta_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(5, [char(plot_loc) '/RoboFly_results/' 'FM_eta_' char(figure_names(i))], 'fig')

    figure(6)
    hFig = figure(6);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 41:60
        hold on
        subplot(3,2,1); plot(time_s(w,:),radtodeg(theta_L_s(w,:)),'Color',[ (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_L [deg]')
    for w = 41:60
        hold on
        subplot(3,2,3); plot(time_s(w,:),90-radtodeg(eta_L_s(w,:)),'Color',[ (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_L [deg]')
    for w = 41:60
        hold on
        subplot(3,2,5); plot(time_s(w,:),-radtodeg(phi_L_s(w,:)),'Color',[ (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_L [deg]')
    for w = 41:60
        hold on
        subplot(3,2,2); plot(time_s(w,:),radtodeg(theta_R_s(w,:)),'Color',[ (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_R [deg]')
    for w = 41:60
        hold on
        subplot(3,2,4); plot(time_s(w,:),90-radtodeg(eta_R_s(w,:)),'Color',[ (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_R [deg]')
    for w = 41:60
        hold on
        subplot(3,2,6); plot(time_s(w,:),-radtodeg(phi_R_s(w,:)),'Color',[ (0.5-(0.5*(w-40))/21) (0.5-(0.5*(w-40))/21) (0.5+(0.5*(w-40))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_R [deg]')
    hold off
    [~,h1] = suplabel(['Wing kinematics ' char(figure_names(i)) ' \eta only'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'Wingkin_eta_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(6, [char(plot_loc) '/RoboFly_results/' 'Wingkin_eta_' char(figure_names(i))], 'fig')

    figure(7)
    hFig = figure(7);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 61:80
        hold on
        subplot(3,2,1); plot(time_s(w,:),FX_s(w,:)-FX_g(w,:),'Color',[ (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fx [N]')
    xlim([5 7])
    ylim([-2 3])
    for w = 61:80
        hold on
        subplot(3,2,3); plot(time_s(w,:),FY_s(w,:)-FY_g(w,:),'Color',[ (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fy [N]')
    xlim([5 7])
    ylim([-1.5 1.5])
    for w = 61:80
        hold on
        subplot(3,2,5); plot(time_s(w,:),-(FZ_s(w,:)-FZ_g(w,:)),'Color',[ (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fz [N]')
    xlim([5 7])
    ylim([-3.5 1])
    for w = 61:80
        hold on
        subplot(3,2,2); plot(time_s(w,:),MX_s(w,:)-MX_g(w,:),'Color',[ (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mx [N*m]')
    xlim([5 7])
    ylim([-0.15 0.25])
    for w = 61:80
        hold on
        subplot(3,2,4); plot(time_s(w,:),MY_s(w,:)-MY_g(w,:),'Color',[ (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('My [N*m]')
    xlim([5 7])
    ylim([-0.25 0.1])
    for w = 61:80
        hold on
        subplot(3,2,6); plot(time_s(w,:),MZ_s(w,:)-MZ_g(w,:),'Color',[ (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mz [N*m]')
    xlim([5 7])
    ylim([-0.25 0.2])
    hold off
    [~,h1] = suplabel(['Forces and Moments ' char(figure_names(i)) ', \phi only'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'FM_phi_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(7, [char(plot_loc) '/RoboFly_results/' 'FM_phi_' char(figure_names(i))], 'fig')

    figure(8)
    hFig = figure(8);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 61:80
        hold on
        subplot(3,2,1); plot(time_s(w,:),radtodeg(theta_L_s(w,:)),'Color',[ (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_L [deg]')
    for w = 61:80
        hold on
        subplot(3,2,3); plot(time_s(w,:),90-radtodeg(eta_L_s(w,:)),'Color',[ (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_L [deg]')
    for w = 61:80
        hold on
        subplot(3,2,5); plot(time_s(w,:),-radtodeg(phi_L_s(w,:)),'Color',[ (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_L [deg]')
    for w = 61:80
        hold on
        subplot(3,2,2); plot(time_s(w,:),radtodeg(theta_R_s(w,:)),'Color',[ (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\theta_R [deg]')
    for w = 61:80
        hold on
        subplot(3,2,4); plot(time_s(w,:),90-radtodeg(eta_R_s(w,:)),'Color',[ (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\eta_R [deg]')
    for w = 61:80
        hold on
        subplot(3,2,6); plot(time_s(w,:),-radtodeg(phi_R_s(w,:)),'Color',[ (0.5-(0.5*(w-60))/21) (0.5-(0.5*(w-60))/21) (0.5+(0.5*(w-60))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('\phi_R [deg]')
    hold off
    [~,h1] = suplabel(['Wing kinematics ' char(figure_names(i)) ' \phi only'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'Wingkin_phi_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(8, [char(plot_loc) '/RoboFly_results/' 'Wingkin_phi_' char(figure_names(i))], 'fig')

    figure(9)
    hFig = figure(9);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    for w = 1:20
        hold on
        subplot(3,2,1); plot(time_s(w,:),FX_g(w,:),'Color',[ (0.5+(0.5*(w))/21) (0.5-(0.5*(w))/21) (0.5-(0.5*(w))/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fx [N]')
    for w = 1:20
        hold on
        subplot(3,2,3); plot(time_s(w,:),FY_g(w,:),'Color',[ (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fy [N]')
    for w = 1:20
        hold on
        subplot(3,2,5); plot(time_s(w,:),-(FZ_g(w,:)),'Color',[ (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Fz [N]')
    for w = 1:20
        hold on
        subplot(3,2,2); plot(time_s(w,:),MX_g(w,:),'Color',[ (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mx [N*m]')
    for w = 1:20
        hold on
        subplot(3,2,4); plot(time_s(w,:),MY_g(w,:),'Color',[ (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) (0.5-(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('My [N*m]')
    for w = 1:20
        hold on
        subplot(3,2,6); plot(time_s(w,:),MZ_g(w,:),'Color',[ (0.5-(0.5*w)/21) (0.5-(0.5*w)/21) (0.5+(0.5*w)/21) ])
        hold off
    end
    xlabel('t/T')
    ylabel('Mz [N*m]')
    hold off
    [~,h1] = suplabel(['Forces and moments ' char(figure_names(i)) ' gravity'], 't');
    set(h1,'FontSize',10)
    print ([char(plot_loc) '/RoboFly_results/' 'FM_gravity_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(9, [char(plot_loc) '/RoboFly_results/' 'FM_gravity_' char(figure_names(i))], 'fig')    
       
    % Find maxima in phi and determine the temporal range of each wingbeat
    % subsequently select the last 4 wingbeats and construct average forces
    % and moments per wingbeat set (left and right).
    
    pks_L_loc = nan(80,7);
    pks_R_loc = nan(80,7);
    
    FX_int = nan(80,3);
    FY_int = nan(80,3);
    FZ_int = nan(80,3);
    MX_int = nan(80,3);
    MY_int = nan(80,3);
    MZ_int = nan(80,3);
    
    FX_g_int = zeros(80,3);
    FY_g_int = zeros(80,3);
    FZ_g_int = zeros(80,3);
    MX_g_int = zeros(80,3);
    MY_g_int = zeros(80,3);
    MZ_g_int = zeros(80,3);

%     FX_int = nan(80,6);
%     FY_int = nan(80,6);
%     FZ_int = nan(80,6);
%     MX_int = nan(80,6);
%     MY_int = nan(80,6);
%     MZ_int = nan(80,6);
%     
%     FX_g_int = zeros(80,6);
%     FY_g_int = zeros(80,6);
%     FZ_g_int = zeros(80,6);
%     MX_g_int = zeros(80,6);
%     MY_g_int = zeros(80,6);
%     MZ_g_int = zeros(80,6);
    
    for n = 1:80
    
        [ ~ , peak_L_loc ] = findpeaks(-phi_L_s(n,:),'minpeakdistance',1000);
        [ ~ , peak_R_loc ] = findpeaks(-phi_R_s(n,:),'minpeakdistance',1000);
        
        pks_L_loc(n,:) = peak_L_loc(1:7);
        pks_R_loc(n,:) = peak_R_loc(1:7);
        

        FX_int(n,:) = [ mean(FX_s(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                        mean(FX_s(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                        mean(FX_s(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5;
                    
        FY_int(n,:) = [ mean(FY_s(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                        mean(FY_s(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                        mean(FY_s(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5;
                    
        FZ_int(n,:) = [ mean(FZ_s(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                        mean(FZ_s(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                        mean(FZ_s(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5;
                    
        MX_int(n,:) = [ mean(MX_s(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                        mean(MX_s(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                        mean(MX_s(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5*(3/230)*1000;
                    
        MY_int(n,:) = [ mean(MY_s(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                        mean(MY_s(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                        mean(MY_s(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5*(3/230)*1000;
                    
        MZ_int(n,:) = [ mean(MZ_s(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                        mean(MZ_s(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                        mean(MZ_s(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5*(3/230)*1000;
                    
        FX_g_int(n,:) = [   mean(FX_g(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                            mean(FX_g(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                            mean(FX_g(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5;
                    
        FY_g_int(n,:) = [   mean(FY_g(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                            mean(FY_g(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                            mean(FY_g(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5;
                    
        FZ_g_int(n,:) = [   mean(FZ_g(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                            mean(FZ_g(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                            mean(FZ_g(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5;
                    
        MX_g_int(n,:) = [   mean(MX_g(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                            mean(MX_g(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                            mean(MX_g(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5*(3/230)*1000;
                    
        MY_g_int(n,:) = [ 	mean(MY_g(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                            mean(MY_g(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                            mean(MY_g(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5*(3/230)*1000;
                    
        MZ_g_int(n,:) = [   mean(MZ_g(n,pks_L_loc(n,4):pks_L_loc(n,5))) ...
                            mean(MZ_g(n,pks_L_loc(n,5):pks_L_loc(n,6))) ...
                            mean(MZ_g(n,pks_L_loc(n,6):pks_L_loc(n,7)))]*2.2429e-5*(3/230)*1000;
       
                
    end

    
%     Fx_back_max         = 0.1864*-0.4015e-4;    % [N]
%     Fx_forward_max      = 0.1864*0.5038e-4;     % [N]
%     Fy_max              = 0.1864*0.4685e-4;     % [N]
%     Fz_up_max           = 0.1864*-0.7458e-4;    % [N]
%     Mx_max              = 0.5508*0.2305e-4;     % [N*mm]
%     My_down_max         = 0.5508*-0.2675e-4;    % [N*mm]
%     My_up_max           = 0.5508*0.2692e-4;     % [N*mm]
%     Mz_max              = 0.5508*0.4288e-4;     % [N*mm]

    Fx_back_max         = 0.1864*-0.4015e-4;    % [N]
    Fx_forward_max      = 0.1864*0.5038e-4;     % [N]
    Fy_max              = 0.1864*0.4685e-4;     % [N]
    Fz_up_max           = 0.1864*-0.6561e-4;    % [N]
    Mx_max              = 0.5508*0.2756e-4;     % [N*mm]
    My_down_max         = 0.5508*-0.3203e-4;    % [N*mm]
    My_up_max           = 0.5508*0.2146e-4;     % [N*mm]
    Mz_max              = 0.5508*0.4255e-4;     % [N*mm]

    FX_0 = ((FX_int(1,1)-FX_g_int(1,1))+(FX_int(1,2)-FX_g_int(1,2))+(FX_int(1,3)-FX_g_int(1,3)))/3;
    FY_0 = ((FY_int(1,1)-FY_g_int(1,1))+(FY_int(1,2)-FY_g_int(1,2))+(FY_int(1,3)-FY_g_int(1,3)))/3;
    FZ_0 = ((FZ_int(1,1)-FZ_g_int(1,1))+(FZ_int(1,2)-FZ_g_int(1,2))+(FZ_int(1,3)-FZ_g_int(1,3)))/3;
    MX_0 = ((MX_int(1,1)-MX_g_int(1,1))+(MX_int(1,2)-MX_g_int(1,2))+(MX_int(1,3)-MX_g_int(1,3)))/3;
    MY_0 = ((MY_int(1,1)-MY_g_int(1,1))+(MY_int(1,2)-MY_g_int(1,2))+(MY_int(1,3)-MY_g_int(1,3)))/3;
    MZ_0 = ((MZ_int(1,1)-MZ_g_int(1,1))+(MZ_int(1,2)-MZ_g_int(1,2))+(MZ_int(1,3)-MZ_g_int(1,3)))/3;
    
    
    figure(10)
    hFig = figure(10);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    if i == 1
        plot(1:20,0:(Fx_back_max/19):(Fx_back_max),'k')
        plot(1:20,FM_qt_ax_back(1,:)-FM_qt_ax_back(1,1),'b')
        plot(1:20,FM_qr_ax_back(1,:)-FM_qr_ax_back(1,1),'m')
    elseif i == 2
        plot(1:20,0:(Fx_forward_max/19):(Fx_forward_max),'k')
        plot(1:20,FM_qt_ax_forward(1,:)-FM_qt_ax_forward(1,1),'b')
        plot(1:20,FM_qr_ax_forward(1,:)-FM_qr_ax_forward(1,1),'m')
    end
    plot(1:20,(((FX_int(1:20,1)-FX_g_int(1:20,1))+(FX_int(1:20,2)-FX_g_int(1:20,2))+(FX_int(1:20,3)-FX_g_int(1:20,3)))/3-FX_0),'r')
    title(char([ 'Fx ' char(figure_names(i))]))
    xlabel('wing kinematic set')
    ylabel('Fx [N]')
    if i == 1 || i == 2
       legend('Predict','QS trans','QS rot','RoboFly')
    else
       legend('Predict','RoboFly')
    end
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict1_Fx_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(10, [char(plot_loc) '/RoboFly_results/' 'predict1_Fx_' char(figure_names(i))], 'fig')
    
    figure(11)
    hFig = figure(11);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    plot(1:20,(((FX_int(1:20,1)-FX_g_int(1:20,1))+(FX_int(1:20,2)-FX_g_int(1:20,2))+(FX_int(1:20,3)-FX_g_int(1:20,3)))/3-FX_0),'k')
    plot(1:20,(((FX_int(21:40,1)-FX_g_int(21:40,1))+(FX_int(21:40,2)-FX_g_int(21:40,2))+(FX_int(21:40,3)-FX_g_int(21:40,3)))/3-FX_0),'r')
    plot(1:20,(((FX_int(41:60,1)-FX_g_int(41:60,1))+(FX_int(41:60,2)-FX_g_int(41:60,2))+(FX_int(41:60,3)-FX_g_int(41:60,3)))/3-FX_0),'g')
    plot(1:20,(((FX_int(61:80,1)-FX_g_int(61:80,1))+(FX_int(61:80,2)-FX_g_int(61:80,2))+(FX_int(61:80,3)-FX_g_int(61:80,3)))/3-FX_0),'b')
    title(char([ 'Fx ' char(figure_names(i))]))
    xlabel('Contributions theta, eta and phi')
    ylabel('Fx [N]')
    legend('full wing kinematics','theta only','eta only','phi only')
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict2_Fx_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(11, [char(plot_loc) '/RoboFly_results/' 'predict2_Fx_' char(figure_names(i))], 'fig')
    
    figure(12)
    hFig = figure(12);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    if i == 3
        plot(1:20,0:(Fy_max/19):(Fy_max),'k')
        plot(1:20,FM_qt_ay(2,:)-FM_qt_ay(2,1),'b')
        plot(1:20,FM_qr_ay(2,:)-FM_qr_ay(2,1),'m')
    end
    plot(1:20,(((FY_int(1:20,1)-FY_g_int(1:20,1))+(FY_int(1:20,2)-FY_g_int(1:20,2))+(FY_int(1:20,3)-FY_g_int(1:20,3)))/3-FY_0),'r')
    title(char([ 'Fy ' char(figure_names(i))]))
    xlabel('wing kinematic set')
    ylabel('Fy [N]')
    if i == 3
       legend('Predict','QS trans','QS rot','RoboFly')
    else
       legend('Predict','RoboFly')
    end
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict1_Fy_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(12, [char(plot_loc) '/RoboFly_results/' 'predict1_Fy_' char(figure_names(i))], 'fig')
    
    figure(13)
    hFig = figure(13);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    plot(1:20,(((FY_int(1:20,1)-FY_g_int(1:20,1))+(FY_int(1:20,2)-FY_g_int(1:20,2))+(FY_int(1:20,3)-FY_g_int(1:20,3)))/3-FY_0),'k')
    plot(1:20,(((FY_int(21:40,1)-FY_g_int(21:40,1))+(FY_int(21:40,2)-FY_g_int(21:40,2))+(FY_int(21:40,3)-FY_g_int(21:40,3)))/3-FY_0),'r')
    plot(1:20,(((FY_int(41:60,1)-FY_g_int(41:60,1))+(FY_int(41:60,2)-FY_g_int(41:60,2))+(FY_int(41:60,3)-FY_g_int(41:60,3)))/3-FY_0),'g')
    plot(1:20,(((FY_int(61:80,1)-FY_g_int(61:80,1))+(FY_int(61:80,2)-FY_g_int(61:80,2))+(FY_int(61:80,3)-FY_g_int(61:80,3)))/3-FY_0),'b')
    title(char([ 'Fy ' char(figure_names(i))]))
    xlabel('Contributions theta, eta and phi')
    ylabel('Fy [N]')
    legend('full wing kinematics','theta only','eta only','phi only')
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict2_Fy_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(13, [char(plot_loc) '/RoboFly_results/' 'predict2_Fy_' char(figure_names(i))], 'fig')
    
    figure(14)
    hFig = figure(14);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    if i == 4
        plot(1:20,0:(Fz_up_max/19):(Fz_up_max),'k')
        plot(1:20,FM_qt_az_up(3,:)-FM_qt_az_up(3,1),'b')
        plot(1:20,FM_qr_az_up(3,:)-FM_qr_az_up(3,1),'m')
    end
    plot(1:20,-(((FZ_int(1:20,1)-FZ_g_int(1:20,1))+(FZ_int(1:20,2)-FZ_g_int(1:20,2))+(FZ_int(1:20,3)-FZ_g_int(1:20,3)))/3-FZ_0),'r')
    title(char([ 'Fz ' char(figure_names(i))]))
    xlabel('wing kinematic set')
    ylabel('Fz [N]')
    if i == 4
       legend('Predict','QS trans','QS rot','RoboFly')
    else
       legend('Predict','RoboFly')
    end
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict1_Fz_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(14, [char(plot_loc) '/RoboFly_results/' 'predict1_Fz_' char(figure_names(i))], 'fig')
    
    figure(15)
    hFig = figure(15);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    plot(1:20,-(((FZ_int(1:20,1)-FZ_g_int(1:20,1))+(FZ_int(1:20,2)-FZ_g_int(1:20,2))+(FZ_int(1:20,3)-FZ_g_int(1:20,3)))/3-FZ_0),'k')
    plot(1:20,-(((FZ_int(21:40,1)-FZ_g_int(21:40,1))+(FZ_int(21:40,2)-FZ_g_int(21:40,2))+(FZ_int(21:40,3)-FZ_g_int(21:40,3)))/3-FZ_0),'r')
    plot(1:20,-(((FZ_int(41:60,1)-FZ_g_int(41:60,1))+(FZ_int(41:60,2)-FZ_g_int(41:60,2))+(FZ_int(41:60,3)-FZ_g_int(41:60,3)))/3-FZ_0),'g')
    plot(1:20,-(((FZ_int(61:80,1)-FZ_g_int(61:80,1))+(FZ_int(61:80,2)-FZ_g_int(61:80,2))+(FZ_int(61:80,3)-FZ_g_int(61:80,3)))/3-FZ_0),'b')
    title(char([ 'Fz ' char(figure_names(i))]))
    xlabel('Contributions theta, eta and phi')
    ylabel('Fz [N]')
    legend('full wing kinematics','theta only','eta only','phi only')
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict2_Fz_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(15, [char(plot_loc) '/RoboFly_results/' 'predict2_Fz_' char(figure_names(i))], 'fig')
       
    figure(16)
    hFig = figure(16);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    if i == 5
        plot(1:20,0:(Mx_max/19):(Mx_max),'k')
        plot(1:20,FM_qt_wx(4,:)-FM_qt_wx(4,1),'b')
        plot(1:20,FM_qr_wx(4,:)-FM_qr_wx(4,1),'m')
    end
    plot(1:20,(((MX_int(1:20,1)-MX_g_int(1:20,1))+(MX_int(1:20,2)-MX_g_int(1:20,2))+(MX_int(1:20,3)-MX_g_int(1:20,3)))/3-MX_0),'r')
    title(char([ 'Mx ' char(figure_names(i))]))
    xlabel('wing kinematic set')
    ylabel('Mx [N*m]')
    if i == 5
       legend('Predict','QS trans','QS rot','RoboFly')
    else
       legend('Predict','RoboFly')
    end
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict1_Mx_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(16, [char(plot_loc) '/RoboFly_results/' 'predict1_Mx_' char(figure_names(i))], 'fig')
    
    figure(17)
    hFig = figure(17);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    plot(1:20,(((MX_int(1:20,1)-MX_g_int(1:20,1))+(MX_int(1:20,2)-MX_g_int(1:20,2))+(MX_int(1:20,3)-MX_g_int(1:20,3)))/3-MX_0),'k')
    plot(1:20,(((MX_int(21:40,1)-MX_g_int(21:40,1))+(MX_int(21:40,2)-MX_g_int(21:40,2))+(MX_int(21:40,3)-MX_g_int(21:40,3)))/3-MX_0),'r')
    plot(1:20,(((MX_int(41:60,1)-MX_g_int(41:60,1))+(MX_int(41:60,2)-MX_g_int(41:60,2))+(MX_int(41:60,3)-MX_g_int(41:60,3)))/3-MX_0),'g')
    plot(1:20,(((MX_int(61:80,1)-MX_g_int(61:80,1))+(MX_int(61:80,2)-MX_g_int(61:80,2))+(MX_int(61:80,3)-MX_g_int(61:80,3)))/3-MX_0),'b')
    title(char([ 'Mx ' char(figure_names(i))]))
    xlabel('Contributions theta, eta and phi')
    ylabel('Mx [N*m]')
    legend('full wing kinematics','theta only','eta only','phi only')
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict2_Mx_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(17, [char(plot_loc) '/RoboFly_results/' 'predict2_Mx_' char(figure_names(i))], 'fig')
    
    figure(18)
    hFig = figure(18);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    if i == 6
        plot(1:20,0:(My_down_max/19):(My_down_max),'k')
        plot(1:20,FM_qt_wy_down(5,:)-FM_qt_wy_down(5,1),'b')
        plot(1:20,FM_qr_wy_down(5,:)-FM_qr_wy_down(5,1),'m')
    end
    if i == 7
        plot(1:20,0:(My_up_max/19):(My_up_max),'k')
        plot(1:20,FM_qt_wy_up(5,:)-FM_qt_wy_up(5,1),'b')
        plot(1:20,FM_qr_wy_up(5,:)-FM_qr_wy_up(5,1),'m')
    end
    plot(1:20,(((MY_int(1:20,1)-MY_g_int(1:20,1))+(MY_int(1:20,2)-MY_g_int(1:20,2))+(MY_int(1:20,3)-MY_g_int(1:20,3)))/3-MY_0),'r')
    title(char([ 'My ' char(figure_names(i))]))
    xlabel('wing kinematic set')
    ylabel('My [N*m]')
    if i == 6 || i == 7
       legend('Predict','QS trans','QS rot','RoboFly')
    else
       legend('Predict','RoboFly')
    end
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict1_My_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(18, [char(plot_loc) '/RoboFly_results/' 'predict1_My_' char(figure_names(i))], 'fig')
    
    figure(19)
    hFig = figure(19);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    plot(1:20,(((MY_int(1:20,1)-MY_g_int(1:20,1))+(MY_int(1:20,2)-MY_g_int(1:20,2))+(MY_int(1:20,3)-MY_g_int(1:20,3)))/3-MY_0),'k')
    plot(1:20,(((MY_int(21:40,1)-MY_g_int(21:40,1))+(MY_int(21:40,2)-MY_g_int(21:40,2))+(MY_int(21:40,3)-MY_g_int(21:40,3)))/3-MY_0),'r')
    plot(1:20,(((MY_int(41:60,1)-MY_g_int(41:60,1))+(MY_int(41:60,2)-MY_g_int(41:60,2))+(MY_int(41:60,3)-MY_g_int(41:60,3)))/3-MY_0),'g')
    plot(1:20,(((MY_int(61:80,1)-MY_g_int(61:80,1))+(MY_int(61:80,2)-MY_g_int(61:80,2))+(MY_int(61:80,3)-MY_g_int(61:80,3)))/3-MY_0),'b')
    title(char([ 'My ' char(figure_names(i))]))
    xlabel('Contributions theta, eta and phi')
    ylabel('My [N*m]')
    legend('full wing kinematics','theta only','eta only','phi only')
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict2_My_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(19, [char(plot_loc) '/RoboFly_results/' 'predict2_My_' char(figure_names(i))], 'fig')
    
    figure(20)
    hFig = figure(20);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    if i == 8
        plot(1:20,0:(Mz_max/19):(Mz_max),'k')
        plot(1:20,FM_qt_wz(6,:)-FM_qt_wz(6,1),'b')
        plot(1:20,FM_qr_wz(6,:)-FM_qr_wz(6,1),'m')
    end
    plot(1:20,-(((MZ_int(1:20,1)-MZ_g_int(1:20,1))+(MZ_int(1:20,2)-MZ_g_int(1:20,2))+(MZ_int(1:20,3)-MZ_g_int(1:20,3)))/3-MZ_0),'r')
    title(char([ 'Mz ' char(figure_names(i))]))
    xlabel('wing kinematic set')
    ylabel('Mz [N*m]')
    if i == 8
       legend('Predict','QS trans','QS rot','RoboFly')
    else
       legend('Predict','RoboFly')
    end
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict1_Mz_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(20, [char(plot_loc) '/RoboFly_results/' 'predict1_Mz_' char(figure_names(i))], 'fig')
    
    figure(21)
    hFig = figure(21);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 600 600]);
    hold on
    plot(1:20,-(((MZ_int(1:20,1)-MZ_g_int(1:20,1))+(MZ_int(1:20,2)-MZ_g_int(1:20,2))+(MZ_int(1:20,3)-MZ_g_int(1:20,3)))/3-MZ_0),'k')
    plot(1:20,-(((MZ_int(21:40,1)-MZ_g_int(21:40,1))+(MZ_int(21:40,2)-MZ_g_int(21:40,2))+(MZ_int(21:40,3)-MZ_g_int(21:40,3)))/3-MZ_0),'r')
    plot(1:20,-(((MZ_int(41:60,1)-MZ_g_int(41:60,1))+(MZ_int(41:60,2)-MZ_g_int(41:60,2))+(MZ_int(41:60,3)-MZ_g_int(41:60,3)))/3-MZ_0),'g')
    plot(1:20,-(((MZ_int(61:80,1)-MZ_g_int(61:80,1))+(MZ_int(61:80,2)-MZ_g_int(61:80,2))+(MZ_int(61:80,3)-MZ_g_int(61:80,3)))/3-MZ_0),'b')
    title(char([ 'MZ ' char(figure_names(i))]))
    xlabel('Contributions theta, eta and phi')
    ylabel('Mz [N*m]')
    legend('full wing kinematics','theta only','eta only','phi only')
    hold off
    print ([char(plot_loc) '/RoboFly_results/' 'predict2_Mz_' char(figure_names(i)) '.eps'] ,'-depsc2');
    saveas(21, [char(plot_loc) '/RoboFly_results/' 'predict2_Mz_' char(figure_names(i))], 'fig')

    close all

end

% % Flowvis experiments:
% 
%     cd('K:/robofly_data_8_7_2013/Johan_06252013/Johan_06272013_FlowVis')
%     
%     temp_dir = dir;
%     
%     dir_names = {temp_dir.name};
% 
%     temp_isdir = [temp_dir.isdir];
%     
%     exp_names1 = [];
%     
%     for j = 1:length(dir_names)
% 
%         if temp_isdir(j) == 0
% 
%            exp_names1 = [exp_names1; dir_names(j)];
% 
%         end
% 
%     end
%     
%     exp_names = {};
%     
%     temp_digits = strfind(char(exp_names1(1)),'_');
%     
%     start_digit = temp_digits(end)+1;
%         
%     for k = 1:length(exp_names1)
%         
%         temp_name = char(exp_names1(k));
%         
%         file_index = str2num(temp_name(start_digit:(end-4)))+1;
%         
%         exp_names(file_index) = exp_names1(k);
%         
%     end
%     
%    
%     
%     Test = [];
% 
%     for m = 1:length(exp_names)
% 
%         temp_file = load((char(exp_names(m))));
% 
%         Test.([ 'test_' int2str(m) ]) = temp_file;
% 
%     end
% 
%     Test_names = fieldnames(Test);
%     
%     nr_vis = length(Test_names);
%     
%     time_vis            = nan(nr_vis,120000);
%     FX_vis              = nan(nr_vis,120000);
%     FY_vis              = nan(nr_vis,120000);
%     FZ_vis              = nan(nr_vis,120000);
%     MX_vis              = nan(nr_vis,120000);
%     MY_vis              = nan(nr_vis,120000);
%     MZ_vis              = nan(nr_vis,120000);
%     theta_L_vis         = nan(nr_vis,120000);
%     eta_L_vis           = nan(nr_vis,120000);
%     phi_L_vis           = nan(nr_vis,120000);
%     theta_R_vis         = nan(nr_vis,120000);
%     eta_R_vis           = nan(nr_vis,120000);
%     phi_R_vis           = nan(nr_vis,120000);
%     freq_vis            = nan(nr_vis,13);
%     wing_length_vis     = nan(nr_vis,13);
% 
%     time_s_vis            = nan(nr_vis,120000);
%     FX_s_vis              = nan(nr_vis,120000);
%     FY_s_vis              = nan(nr_vis,120000);
%     FZ_s_vis              = nan(nr_vis,120000);
%     MX_s_vis              = nan(nr_vis,120000);
%     MY_s_vis              = nan(nr_vis,120000);
%     MZ_s_vis              = nan(nr_vis,120000);
%     
%     for n = 1:nr_vis
%         
%         test_end = length(Test.(char(Test_names(n))).t);
%         
%         time_vis(n,1:test_end)            = Test.(char(Test_names(n))).t;
%         FX_vis(n,1:test_end)              = Test.(char(Test_names(n))).ft(:,1);
%         FY_vis(n,1:test_end)              = Test.(char(Test_names(n))).ft(:,3);
%         FZ_vis(n,1:test_end)              = Test.(char(Test_names(n))).ft(:,2);
%         MX_vis(n,1:test_end)              = Test.(char(Test_names(n))).ft(:,4);
%         MY_vis(n,1:test_end)              = Test.(char(Test_names(n))).ft(:,6);
%         MZ_vis(n,1:test_end)              = Test.(char(Test_names(n))).ft(:,5);
%         theta_L_vis(n,1:test_end)         = Test.(char(Test_names(n))).kine(:,3);
%         eta_L_vis(n,1:test_end)           = Test.(char(Test_names(n))).kine(:,1);
%         phi_L_vis(n,1:test_end)           = Test.(char(Test_names(n))).kine(:,5);
%         theta_R_vis(n,1:test_end)         = Test.(char(Test_names(n))).kine(:,2);
%         eta_R_vis(n,1:test_end)           = Test.(char(Test_names(n))).kine(:,4);
%         phi_R_vis(n,1:test_end)           = Test.(char(Test_names(n))).kine(:,6);
%         freq_vis(n,:)              = Test.(char(Test_names(n))).freq;
%         wing_length_vis(n,:)       = Test.(char(Test_names(n))).wing_length;
%         
%         time_s_vis(n,1:test_end)            = time_vis(n,1:test_end)./((1/13)*time_vis(n,test_end));
%         FX_s_vis(n,1:test_end)              = 2.613e-5*flipud(filter(ones(1,50)/50,1,flipud(filter(ones(1,50)/50,1,FX_vis(n,1:test_end)))));
%         FY_s_vis(n,1:test_end)              = 2.613e-5*flipud(filter(ones(1,50)/50,1,flipud(filter(ones(1,50)/50,1,FY_vis(n,1:test_end)))));
%         FZ_s_vis(n,1:test_end)              = 2.613e-5*flipud(filter(ones(1,50)/50,1,flipud(filter(ones(1,50)/50,1,FZ_vis(n,1:test_end)))));
%         MX_s_vis(n,1:test_end)              = 1e-3*2.613e-5*flipud(filter(ones(1,50)/50,1,flipud(filter(ones(1,50)/50,1,MX_vis(n,1:test_end)))));
%         MY_s_vis(n,1:test_end)              = 1e-3*2.613e-5*flipud(filter(ones(1,50)/50,1,flipud(filter(ones(1,50)/50,1,MY_vis(n,1:test_end)))));
%         MZ_s_vis(n,1:test_end)              = 1e-3*2.613e-5*flipud(filter(ones(1,50)/50,1,flipud(filter(ones(1,50)/50,1,MZ_vis(n,1:test_end)))));
%         
%         
%     end
%     
%     figure()
%     hold on
%     for n = 1:nr_vis 
%     plot(time_s_vis(n,:),FX_s_vis(n,:),'Color',[ (0.5+(0.5*n)/6) (0.5-(0.5*n)/6) (0.5-(0.5*n)/6) ])
%     end
%     hold off
%     
%     figure()
%     hold on
%     for n = 1:nr_vis 
%     plot(time_s_vis(n,:),FY_s_vis(n,:),'Color',[ (0.5-(0.5*n)/6) (0.5+(0.5*n)/6) (0.5-(0.5*n)/6) ])
%     end
%     hold off
%     
%     figure()
%     hold on
%     for n = 1:nr_vis 
%     plot(time_s_vis(n,:),FZ_s_vis(n,:),'Color',[ (0.5-(0.5*n)/6) (0.5-(0.5*n)/6) (0.5+(0.5*n)/6) ])
%     end
%     hold off
%         
%     figure()
%     hold on
%     for n = 1:nr_vis 
%     plot(time_s_vis(n,:),MX_s_vis(n,:),'Color',[ (0.5+(0.5*n)/6) (0.5-(0.5*n)/6) (0.5-(0.5*n)/6) ])
%     end
%     hold off
%         
%     figure()
%     hold on
%     for n = 1:nr_vis 
%     plot(time_s_vis(n,:),MY_s_vis(n,:),'Color',[ (0.5-(0.5*n)/6) (0.5+(0.5*n)/6) (0.5-(0.5*n)/6) ])
%     end
%     hold off
%     
%     figure()
%     hold on
%     for n = 1:nr_vis 
%     plot(time_s_vis(n,:),MZ_s_vis(n,:),'Color',[ (0.5-(0.5*n)/6) (0.5-(0.5*n)/6) (0.5+(0.5*n)/6) ])
%     end
%     hold off
    