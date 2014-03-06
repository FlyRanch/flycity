clear all;
close all;
clc;

% Program to analyze Robofly data:


add_paths = {'W:/RoboFly_data_30_05_2013'};


addpath(char(add_paths(1)))

folder_names = [];

temp_dir = dir;

dir_names = {temp_dir.name};

temp_isdir = [temp_dir.isdir];

for i = 1:length(dir_names)
           
    if temp_isdir(i) == 0
       
        folder_names = [folder_names; dir_names(i)];
        
    end
    
end

Test = [];

for i = 1:length(folder_names)
    
    temp_file = load((char(folder_names(i))));
    
    temp_name = char(folder_names(i));
    
    pos1 = findstr(temp_name,'_');
    
    pos2 = findstr(temp_name,'.');
    
    temp_nr = str2num(temp_name((pos1(end)+1):(pos2-1)));
    
    Test.([char(temp_file.exp_type) '_' int2str(temp_nr+1) ]) = temp_file;
    
end

Test_names = fieldnames(Test);



Wingbeats = {};

Wingbeats.wb_glob_1.t = Test.glob_test_1.t;

Wingbeats.wb_glob_1.theta_L = Test.glob_test_1.kine(:,3);
Wingbeats.wb_glob_1.eta_L = Test.glob_test_1.kine(:,1);
Wingbeats.wb_glob_1.phi_L = Test.glob_test_1.kine(:,5);

Wingbeats.wb_glob_1.theta_R = Test.glob_test_1.kine(:,2);
Wingbeats.wb_glob_1.eta_R = Test.glob_test_1.kine(:,4);
Wingbeats.wb_glob_1.phi_R = Test.glob_test_1.kine(:,6);

Wingbeats.wb_glob_1.Fx = Test.glob_test_1.ft(:,1);
Wingbeats.wb_glob_1.Fy = Test.glob_test_1.ft(:,3);
Wingbeats.wb_glob_1.Fz = Test.glob_test_1.ft(:,2);

Wingbeats.wb_glob_1.Mx = Test.glob_test_1.ft(:,4);
Wingbeats.wb_glob_1.My = Test.glob_test_1.ft(:,6);
Wingbeats.wb_glob_1.Mz = Test.glob_test_1.ft(:,5);

Wingbeats.wb_glob_1.f_robo = Test.glob_test_1.freq;

Wingbeats.wb_glob_1.wing_l = Test.glob_test_1.wing_length;



Wingbeats.wb_glob_2.t = Test.glob_test_2.t;

Wingbeats.wb_glob_2.theta_L = Test.glob_test_2.kine(:,3);
Wingbeats.wb_glob_2.eta_L = Test.glob_test_2.kine(:,1);
Wingbeats.wb_glob_2.phi_L = Test.glob_test_2.kine(:,5);

Wingbeats.wb_glob_2.theta_R = Test.glob_test_2.kine(:,2);
Wingbeats.wb_glob_2.eta_R = Test.glob_test_2.kine(:,4);
Wingbeats.wb_glob_2.phi_R = Test.glob_test_2.kine(:,6);

Wingbeats.wb_glob_2.Fx = Test.glob_test_2.ft(:,1);
Wingbeats.wb_glob_2.Fy = Test.glob_test_2.ft(:,3);
Wingbeats.wb_glob_2.Fz = Test.glob_test_2.ft(:,2);

Wingbeats.wb_glob_2.Mx = Test.glob_test_2.ft(:,4);
Wingbeats.wb_glob_2.My = Test.glob_test_2.ft(:,6);
Wingbeats.wb_glob_2.Mz = Test.glob_test_2.ft(:,5);

Wingbeats.wb_glob_2.f_robo = Test.glob_test_2.freq;

Wingbeats.wb_glob_2.wing_l = Test.glob_test_2.wing_length;



Wingbeats.wb_glob_3.t = Test.glob_test_3.t;

Wingbeats.wb_glob_3.theta_L = Test.glob_test_3.kine(:,3);
Wingbeats.wb_glob_3.eta_L = Test.glob_test_3.kine(:,1);
Wingbeats.wb_glob_3.phi_L = Test.glob_test_3.kine(:,5);

Wingbeats.wb_glob_3.theta_R = Test.glob_test_3.kine(:,2);
Wingbeats.wb_glob_3.eta_R = Test.glob_test_3.kine(:,4);
Wingbeats.wb_glob_3.phi_R = Test.glob_test_3.kine(:,6);

Wingbeats.wb_glob_3.Fx = Test.glob_test_3.ft(:,1);
Wingbeats.wb_glob_3.Fy = Test.glob_test_3.ft(:,3);
Wingbeats.wb_glob_3.Fz = Test.glob_test_3.ft(:,2);

Wingbeats.wb_glob_3.Mx = Test.glob_test_3.ft(:,4);
Wingbeats.wb_glob_3.My = Test.glob_test_3.ft(:,6);
Wingbeats.wb_glob_3.Mz = Test.glob_test_3.ft(:,5);

Wingbeats.wb_glob_3.f_robo = Test.glob_test_3.freq;

Wingbeats.wb_glob_3.wing_l = Test.glob_test_3.wing_length;




Wingbeats.wb_sym_1.t = Test.sym_test_1.t;

Wingbeats.wb_sym_1.theta_L = Test.sym_test_1.kine(:,3);
Wingbeats.wb_sym_1.eta_L = Test.sym_test_1.kine(:,1);
Wingbeats.wb_sym_1.phi_L = Test.sym_test_1.kine(:,5);

Wingbeats.wb_sym_1.theta_R = Test.sym_test_1.kine(:,2);
Wingbeats.wb_sym_1.eta_R = Test.sym_test_1.kine(:,4);
Wingbeats.wb_sym_1.phi_R = Test.sym_test_1.kine(:,6);

Wingbeats.wb_sym_1.Fx = Test.sym_test_1.ft(:,1);
Wingbeats.wb_sym_1.Fy = Test.sym_test_1.ft(:,3);
Wingbeats.wb_sym_1.Fz = Test.sym_test_1.ft(:,2);

Wingbeats.wb_sym_1.Mx = Test.sym_test_1.ft(:,4);
Wingbeats.wb_sym_1.My = Test.sym_test_1.ft(:,6);
Wingbeats.wb_sym_1.Mz = Test.sym_test_1.ft(:,5);

Wingbeats.wb_sym_1.f_robo = Test.sym_test_1.freq;

Wingbeats.wb_sym_1.wing_l = Test.sym_test_1.wing_length;




Wingbeats.wb_sym_2.t = Test.sym_test_2.t;

Wingbeats.wb_sym_2.theta_L = Test.sym_test_2.kine(:,3);
Wingbeats.wb_sym_2.eta_L = Test.sym_test_2.kine(:,1);
Wingbeats.wb_sym_2.phi_L = Test.sym_test_2.kine(:,5);

Wingbeats.wb_sym_2.theta_R = Test.sym_test_2.kine(:,2);
Wingbeats.wb_sym_2.eta_R = Test.sym_test_2.kine(:,4);
Wingbeats.wb_sym_2.phi_R = Test.sym_test_2.kine(:,6);

Wingbeats.wb_sym_2.Fx = Test.sym_test_2.ft(:,1);
Wingbeats.wb_sym_2.Fy = Test.sym_test_2.ft(:,3);
Wingbeats.wb_sym_2.Fz = Test.sym_test_2.ft(:,2);

Wingbeats.wb_sym_2.Mx = Test.sym_test_2.ft(:,4);
Wingbeats.wb_sym_2.My = Test.sym_test_2.ft(:,6);
Wingbeats.wb_sym_2.Mz = Test.sym_test_2.ft(:,5);

Wingbeats.wb_sym_2.f_robo = Test.sym_test_2.freq;

Wingbeats.wb_sym_2.wing_l = Test.sym_test_2.wing_length;




Wingbeats.wb_sym_3.t = Test.sym_test_3.t;

Wingbeats.wb_sym_3.theta_L = Test.sym_test_3.kine(:,3);
Wingbeats.wb_sym_3.eta_L = Test.sym_test_3.kine(:,1);
Wingbeats.wb_sym_3.phi_L = Test.sym_test_3.kine(:,5);

Wingbeats.wb_sym_3.theta_R = Test.sym_test_3.kine(:,2);
Wingbeats.wb_sym_3.eta_R = Test.sym_test_3.kine(:,4);
Wingbeats.wb_sym_3.phi_R = Test.sym_test_3.kine(:,6);

Wingbeats.wb_sym_3.Fx = Test.sym_test_3.ft(:,1);
Wingbeats.wb_sym_3.Fy = Test.sym_test_3.ft(:,3);
Wingbeats.wb_sym_3.Fz = Test.sym_test_3.ft(:,2);

Wingbeats.wb_sym_3.Mx = Test.sym_test_3.ft(:,4);
Wingbeats.wb_sym_3.My = Test.sym_test_3.ft(:,6);
Wingbeats.wb_sym_3.Mz = Test.sym_test_3.ft(:,5);

Wingbeats.wb_sym_3.f_robo = Test.sym_test_3.freq;

Wingbeats.wb_sym_3.wing_l = Test.sym_test_3.wing_length;




Wingbeats.wb_sym_4.t = Test.sym_test_4.t;

Wingbeats.wb_sym_4.theta_L = Test.sym_test_4.kine(:,3);
Wingbeats.wb_sym_4.eta_L = Test.sym_test_4.kine(:,1);
Wingbeats.wb_sym_4.phi_L = Test.sym_test_4.kine(:,5);

Wingbeats.wb_sym_4.theta_R = Test.sym_test_4.kine(:,2);
Wingbeats.wb_sym_4.eta_R = Test.sym_test_4.kine(:,4);
Wingbeats.wb_sym_4.phi_R = Test.sym_test_4.kine(:,6);

Wingbeats.wb_sym_4.Fx = Test.sym_test_4.ft(:,1);
Wingbeats.wb_sym_4.Fy = Test.sym_test_4.ft(:,3);
Wingbeats.wb_sym_4.Fz = Test.sym_test_4.ft(:,2);

Wingbeats.wb_sym_4.Mx = Test.sym_test_4.ft(:,4);
Wingbeats.wb_sym_4.My = Test.sym_test_4.ft(:,6);
Wingbeats.wb_sym_4.Mz = Test.sym_test_4.ft(:,5);

Wingbeats.wb_sym_4.f_robo = Test.sym_test_4.freq;

Wingbeats.wb_sym_4.wing_l = Test.sym_test_4.wing_length;




Wingbeats.wb_sym_5.t = Test.sym_test_5.t;

Wingbeats.wb_sym_5.theta_L = Test.sym_test_5.kine(:,3);
Wingbeats.wb_sym_5.eta_L = Test.sym_test_5.kine(:,1);
Wingbeats.wb_sym_5.phi_L = Test.sym_test_5.kine(:,5);

Wingbeats.wb_sym_5.theta_R = Test.sym_test_5.kine(:,2);
Wingbeats.wb_sym_5.eta_R = Test.sym_test_5.kine(:,4);
Wingbeats.wb_sym_5.phi_R = Test.sym_test_5.kine(:,6);

Wingbeats.wb_sym_5.Fx = Test.sym_test_5.ft(:,1);
Wingbeats.wb_sym_5.Fy = Test.sym_test_5.ft(:,3);
Wingbeats.wb_sym_5.Fz = Test.sym_test_5.ft(:,2);

Wingbeats.wb_sym_5.Mx = Test.sym_test_5.ft(:,4);
Wingbeats.wb_sym_5.My = Test.sym_test_5.ft(:,6);
Wingbeats.wb_sym_5.Mz = Test.sym_test_5.ft(:,5);

Wingbeats.wb_sym_5.f_robo = Test.sym_test_5.freq;

Wingbeats.wb_sym_5.wing_l = Test.sym_test_5.wing_length;




Wingbeats.wb_asym_1.t = Test.asym_test1_1.t;

Wingbeats.wb_asym_1.theta_L = Test.asym_test1_1.kine(:,3);
Wingbeats.wb_asym_1.eta_L = Test.asym_test1_1.kine(:,1);
Wingbeats.wb_asym_1.phi_L = Test.asym_test1_1.kine(:,5);

Wingbeats.wb_asym_1.theta_R = Test.asym_test1_1.kine(:,2);
Wingbeats.wb_asym_1.eta_R = Test.asym_test1_1.kine(:,4);
Wingbeats.wb_asym_1.phi_R = Test.asym_test1_1.kine(:,6);

Wingbeats.wb_asym_1.Fx = Test.asym_test1_1.ft(:,1);
Wingbeats.wb_asym_1.Fy = Test.asym_test1_1.ft(:,3);
Wingbeats.wb_asym_1.Fz = Test.asym_test1_1.ft(:,2);

Wingbeats.wb_asym_1.Mx = Test.asym_test1_1.ft(:,4);
Wingbeats.wb_asym_1.My = Test.asym_test1_1.ft(:,6);
Wingbeats.wb_asym_1.Mz = Test.asym_test1_1.ft(:,5);

Wingbeats.wb_asym_1.f_robo = Test.asym_test1_1.freq;

Wingbeats.wb_asym_1.wing_l = Test.asym_test1_1.wing_length;



Wingbeats.wb_asym_2.t = Test.asym_test1_2.t;

Wingbeats.wb_asym_2.theta_L = Test.asym_test1_2.kine(:,3);
Wingbeats.wb_asym_2.eta_L = Test.asym_test1_2.kine(:,1);
Wingbeats.wb_asym_2.phi_L = Test.asym_test1_2.kine(:,5);

Wingbeats.wb_asym_2.theta_R = Test.asym_test1_2.kine(:,2);
Wingbeats.wb_asym_2.eta_R = Test.asym_test1_2.kine(:,4);
Wingbeats.wb_asym_2.phi_R = Test.asym_test1_2.kine(:,6);

Wingbeats.wb_asym_2.Fx = Test.asym_test1_2.ft(:,1);
Wingbeats.wb_asym_2.Fy = Test.asym_test1_2.ft(:,3);
Wingbeats.wb_asym_2.Fz = Test.asym_test1_2.ft(:,2);

Wingbeats.wb_asym_2.Mx = Test.asym_test1_2.ft(:,4);
Wingbeats.wb_asym_2.My = Test.asym_test1_2.ft(:,6);
Wingbeats.wb_asym_2.Mz = Test.asym_test1_2.ft(:,5);

Wingbeats.wb_asym_2.f_robo = Test.asym_test1_2.freq;

Wingbeats.wb_asym_2.wing_l = Test.asym_test1_2.wing_length;




Wingbeats.wb_asym_3.t = Test.asym_test1_3.t;

Wingbeats.wb_asym_3.theta_L = Test.asym_test1_3.kine(:,3);
Wingbeats.wb_asym_3.eta_L = Test.asym_test1_3.kine(:,1);
Wingbeats.wb_asym_3.phi_L = Test.asym_test1_3.kine(:,5);

Wingbeats.wb_asym_3.theta_R = Test.asym_test1_3.kine(:,2);
Wingbeats.wb_asym_3.eta_R = Test.asym_test1_3.kine(:,4);
Wingbeats.wb_asym_3.phi_R = Test.asym_test1_3.kine(:,6);

Wingbeats.wb_asym_3.Fx = Test.asym_test1_3.ft(:,1);
Wingbeats.wb_asym_3.Fy = Test.asym_test1_3.ft(:,3);
Wingbeats.wb_asym_3.Fz = Test.asym_test1_3.ft(:,2);

Wingbeats.wb_asym_3.Mx = Test.asym_test1_3.ft(:,4);
Wingbeats.wb_asym_3.My = Test.asym_test1_3.ft(:,6);
Wingbeats.wb_asym_3.Mz = Test.asym_test1_3.ft(:,5);

Wingbeats.wb_asym_3.f_robo = Test.asym_test1_3.freq;

Wingbeats.wb_asym_3.wing_l = Test.asym_test1_3.wing_length;




Wingbeats.wb_asym_4.t = Test.asym_test1_4.t;

Wingbeats.wb_asym_4.theta_L = Test.asym_test1_4.kine(:,3);
Wingbeats.wb_asym_4.eta_L = Test.asym_test1_4.kine(:,1);
Wingbeats.wb_asym_4.phi_L = Test.asym_test1_4.kine(:,5);

Wingbeats.wb_asym_4.theta_R = Test.asym_test1_4.kine(:,2);
Wingbeats.wb_asym_4.eta_R = Test.asym_test1_4.kine(:,4);
Wingbeats.wb_asym_4.phi_R = Test.asym_test1_4.kine(:,6);

Wingbeats.wb_asym_4.Fx = Test.asym_test1_4.ft(:,1);
Wingbeats.wb_asym_4.Fy = Test.asym_test1_4.ft(:,3);
Wingbeats.wb_asym_4.Fz = Test.asym_test1_4.ft(:,2);

Wingbeats.wb_asym_4.Mx = Test.asym_test1_4.ft(:,4);
Wingbeats.wb_asym_4.My = Test.asym_test1_4.ft(:,6);
Wingbeats.wb_asym_4.Mz = Test.asym_test1_4.ft(:,5);

Wingbeats.wb_asym_4.f_robo = Test.asym_test1_4.freq;

Wingbeats.wb_asym_4.wing_l = Test.asym_test1_4.wing_length;




Wingbeats.wb_asym_5.t = Test.asym_test1_5.t;

Wingbeats.wb_asym_5.theta_L = Test.asym_test1_5.kine(:,3);
Wingbeats.wb_asym_5.eta_L = Test.asym_test1_5.kine(:,1);
Wingbeats.wb_asym_5.phi_L = Test.asym_test1_5.kine(:,5);

Wingbeats.wb_asym_5.theta_R = Test.asym_test1_5.kine(:,2);
Wingbeats.wb_asym_5.eta_R = Test.asym_test1_5.kine(:,4);
Wingbeats.wb_asym_5.phi_R = Test.asym_test1_5.kine(:,6);

Wingbeats.wb_asym_5.Fx = Test.asym_test1_5.ft(:,1);
Wingbeats.wb_asym_5.Fy = Test.asym_test1_5.ft(:,3);
Wingbeats.wb_asym_5.Fz = Test.asym_test1_5.ft(:,2);

Wingbeats.wb_asym_5.Mx = Test.asym_test1_5.ft(:,4);
Wingbeats.wb_asym_5.My = Test.asym_test1_5.ft(:,6);
Wingbeats.wb_asym_5.Mz = Test.asym_test1_5.ft(:,5);

Wingbeats.wb_asym_5.f_robo = Test.asym_test1_5.freq;

Wingbeats.wb_asym_5.wing_l = Test.asym_test1_5.wing_length;




Wingbeats_names = fieldnames(Wingbeats);


for i = 1:length(Wingbeats_names)
    
    temp_wb = Wingbeats.(char(Wingbeats_names(i)));
    
    meas_start = round(length(temp_wb.t)*(4/7));
    meas_end = length(temp_wb.t);
    meas_r = meas_start:meas_end;
    
    theta_L = radtodeg(temp_wb.theta_L(meas_r));
    eta_L = radtodeg(temp_wb.eta_L(meas_r));
    phi_L = radtodeg(temp_wb.phi_L(meas_r));
    
    theta_R = radtodeg(temp_wb.theta_R(meas_r));
    eta_R = radtodeg(temp_wb.eta_R(meas_r));
    phi_R = radtodeg(temp_wb.phi_R(meas_r));
    
%     t = temp_wb.t(meas_r)*(115/16)*(temp_wb.wing_l/0.23)^2;

    t = meas_r./((1/3)*length(meas_r));
    
    Fx = 2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb.Fx(meas_r)))));
    Fy = 2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb.Fy(meas_r)))));
    Fz = 2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb.Fz(meas_r)))));
    
    Mx = 1e-3*2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb.Mx(meas_r)))));
    My = 1e-3*2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb.My(meas_r)))));
    Mz = 1e-3*2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb.Mz(meas_r)))));
    
    
    temp_wb_glob = Wingbeats.(char(Wingbeats_names(1)));

    meas_start_glob = round(length(temp_wb_glob.t)*(4/7));
    meas_end_glob = length(temp_wb_glob.t);
    meas_r_glob = meas_start_glob:meas_end_glob;

    theta_L_glob = radtodeg(temp_wb_glob.theta_L(meas_r_glob));
    eta_L_glob = radtodeg(temp_wb_glob.eta_L(meas_r_glob));
    phi_L_glob = radtodeg(temp_wb_glob.phi_L(meas_r_glob));

    theta_R_glob = radtodeg(temp_wb_glob.theta_R(meas_r_glob));
    eta_R_glob = radtodeg(temp_wb_glob.eta_R(meas_r_glob));
    phi_R_glob = radtodeg(temp_wb_glob.phi_R(meas_r_glob));

%     t_glob = temp_wb_glob.t(meas_r_glob)*(115/16)*(temp_wb_glob.wing_l/0.23)^2;

    t_glob = meas_r_glob./((1/3)*length(meas_r_glob));

    Fx_glob = 2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb_glob.Fx(meas_r_glob)))));
    Fy_glob = 2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb_glob.Fy(meas_r_glob)))));
    Fz_glob = 2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb_glob.Fz(meas_r_glob)))));

    Mx_glob = 1e-3*2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb_glob.Mx(meas_r_glob)))));
    My_glob = 1e-3*2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb_glob.My(meas_r_glob)))));
    Mz_glob = 1e-3*2.613e-5*flipud(filter(ones(1,150)/150,1,flipud(filter(ones(1,150)/150,1,temp_wb_glob.Mz(meas_r_glob)))));

        
%         figure()
%         hold on
%         subplot(2,2,1); plot(t,theta_L,t,eta_L,t,phi_L,t_glob,theta_L_glob,t_glob,eta_L_glob,t_glob,phi_L_glob)
%         title('Left wing kinematics')
%         xlabel('time [s]')
%         ylabel('Angle [deg]')
%         legend('\theta','\eta','\phi','\theta_{m}','\eta_{m}','\phi_{m}')
%         subplot(2,2,3); plot(t,Fx,t,Fy,t,Fz,t_glob,Fx_glob,t_glob,Fy_glob,t_glob,Fz_glob)
%         title('Left forces')
%         xlabel('time [s]')
%         ylabel('Force [N]')
%         legend('F_x','F_y','F_z','F_{xm}','F_{ym}','F_{zm}')
%         subplot(2,2,4); plot(t,Mx,t,My,t,Mz,t_glob,Mx_glob,t_glob,My_glob,t_glob,Mz_glob)
%         title('Left moments')
%         xlabel('time [s]')
%         ylabel('Moment [N*mm]')
%         legend('M_x','M_y','M_z','M_{xm}','M_{ym}','M_{zm}')
%         subplot(2,2,2); plot(t,theta_R,t,eta_R,t,phi_R,t_glob,theta_R_glob,t_glob,eta_R_glob,t_glob,phi_R_glob)
%         title('Right wing kinematics')
%         xlabel('time [s]')
%         ylabel('Angle [deg]')
%         legend('\theta','\eta','\phi','\theta_{m}','\eta_{m}','\phi_{m}')
%         hold off
        
        figure()
        hold on
        subplot(3,3,1); hold on
        plot(t,theta_L,'r',t,theta_R,'b')
        plot(t_glob,theta_L_glob,'m')
        hold off
        title('\theta')
        xlabel('wingbeats')
        ylabel('Angle [deg]')
        legend('left','right','mean')
        subplot(3,3,2); hold on
        plot(t,eta_L,'r',t,eta_R,'b')
        plot(t_glob,eta_L_glob,'m')
        hold off
        title('\eta')
        xlabel('wingbeats')
        ylabel('Angle [deg]')
        legend('left','right','mean')
        subplot(3,3,3); hold on
        plot(t,phi_L,'r',t,phi_R,'b')
        plot(t_glob,phi_L_glob,'m')
        hold off
        title('\phi')
        xlabel('wingbeats')
        ylabel('Angle [deg]')
        legend('left','right','mean')
        subplot(3,3,4); hold on
        plot(t,Fx,'r')
        plot(t_glob,Fx_glob,'Color',[0.5 0.5 0.5])
        hold off
        title('F_x')
        xlabel('wingbeats')
        ylabel('Force [N]')
        legend('current','mean')
        subplot(3,3,5); hold on
        plot(t,Fy,'r')
        plot(t_glob,Fy_glob,'Color',[0.5 0.5 0.5])
        hold off
        title('F_y')
        xlabel('wingbeats')
        ylabel('Force [N]')
        legend('current','mean')
        subplot(3,3,6); hold on
        plot(t,Fz,'r')
        plot(t_glob,Fz_glob,'Color',[0.5 0.5 0.5])
        hold off
        title('F_z')
        xlabel('wingbeats')
        ylabel('Force [N]')
        legend('current','mean')
        subplot(3,3,7); hold on
        plot(t,Mx,'r')
        plot(t_glob,Mx_glob,'Color',[0.5 0.5 0.5])
        hold off
        title('M_x')
        xlabel('wingbeats')
        ylabel('Moment [N*mm]')
        legend('current','mean')
        subplot(3,3,8); hold on
        plot(t,My,'r')
        plot(t_glob,My_glob,'Color',[0.5 0.5 0.5])
        hold off
        title('M_y')
        xlabel('wingbeats')
        ylabel('Moment [N*mm]')
        legend('current','mean')
        subplot(3,3,9); hold on
        plot(t,Mz,'r')
        plot(t_glob,Mz_glob,'Color',[0.5 0.5 0.5])
        hold off
        title('M_z')
        xlabel('wingbeats')
        ylabel('Moment [N*mm]')
        legend('current','mean')
        hold off

%     end
    
end