function [temp] = pathDB3_nosave(settings,pathDB)
% addpath('/home/florian/Dropbox/WORK/flytracker/flytracker')
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker')
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/mex/')
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/core/')
addpath('/home/florian/Dropbox/FlyCity/WORK/flytracker/results/')

addpath('/home/matt/Dropbox/WORK/flytracker')
addpath('/home/matt/Dropbox/WORK/flytracker/mex/')
addpath('/home/matt/Dropbox/WORK/flytracker/core/')
addpath('/home/matt/Dropbox/WORK/flytracker/results/')

savefile = 'pathDB3.mat';

%Function that creates file pathDB3.mat with filtered Q_body correction for
%Q_wings and angular rotations in the body frame of reference.

qL1_filt1 = nan(size(pathDB.x));
qL2_filt1 = nan(size(pathDB.x));
qL3_filt1 = nan(size(pathDB.x));
qL4_filt1 = nan(size(pathDB.x));

qR1_filt1 = nan(size(pathDB.x));
qR2_filt1 = nan(size(pathDB.x));
qR3_filt1 = nan(size(pathDB.x));
qR4_filt1 = nan(size(pathDB.x));

alfa_body = nan(size(pathDB.x));
beta_body = nan(size(pathDB.x));

roll_body = nan(size(pathDB.x));
pitch_body = nan(size(pathDB.x));
yaw_body = nan(size(pathDB.x));

u_body = nan(size(pathDB.x));
v_body = nan(size(pathDB.x));
w_body = nan(size(pathDB.x));

ax_body = nan(size(pathDB.x));
ay_body = nan(size(pathDB.x));
az_body = nan(size(pathDB.x));

body_length = nan(1,size(pathDB.x,2));
wing_length = nan(1,size(pathDB.x,2));
joint_pos_L = nan(3,size(pathDB.x,2));
joint_pos_R = nan(3,size(pathDB.x,2));


for i=1:size(pathDB.qb1_filt,2)
    
    % start and stop point for the measurements
    start = find(isnan(pathDB.qb1(:,i))==0, 1 );
    stop = find(isnan(pathDB.qb1(:,i))==0, 1, 'last' );
    
%     addpath(char(settings.path_names(7)));
%     addpath(char(settings.path_names(8)));
    
    %Calculate the body length and wing lengths
        
%         alldirs=dir;
        alldirs=settings.sequence_names;
        
        %Calculate the postions of the left and right wing hinges w.r.t. the bodyframe
        load_name = [alldirs{i}, '/flytracks/ManualFit_flytracks'];
%         load_name = strcat(alldirs(i+2).name, '/flytracks/ManualFit_flytracks');
        load(load_name);

        BL = ManualFit.params.bodyscale*(ManualFit.params.bodylen+ManualFit.params.headlen);
        
        params.wing = ManualFit.params.wing;
        params.winglen = ManualFit.params.winglen;
        
        % Assign model parameters
        PAR.params = ManualFit.params;
        PAR.DLT = ManualFit.DLT;
        PAR.cam = ManualFit.cam;
        
       
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
        
        clear ManualFit
        
%         wing_l = (abs(max(ManualFit.params.wing(:,1)))+abs(min(ManualFit.params.wing(:,1))))*ManualFit.params.wingscale;

        SOLN = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1];

        clear flymodQ
        [x,y,z] = flymodQ(SOLN,PAR.params,PAR);
        for j = 1:length(x);
            PAR.modsample(j) = size(x{j},1);
        end
        
        ywL = y{2};
        ywR = y{3};

 
        wing_l_L = max(ywL(:,1))-min(ywL(:,1));
        wing_l_R = max(ywR(:,1))-min(ywR(:,1));
        
        if wing_l_L == wing_l_R
            wing_l = wing_l_L;
        else
            'error: non equal winglengths'
        end
        
        body_length(i) = BL;
        wing_length(i) = wing_l;
        joint_pos_L(:,i) = BL.*([0.2021; -0.1055; -0.1477]);
        joint_pos_R(:,i) = BL.*([0.2021; 0.1055; -0.1477]);
%         
%         figure()
%         plot3(x{1},y{1},z{1},'b')
%         hold on
%         plot3(x{2},y{2},z{2},'r')
%         plot3(x{3},y{3},z{3},'g')
%         plot3([joint_pos_L(1,i) joint_pos_L(1,i)], [joint_pos_L(2,i) joint_pos_L(2,i)-wing_l], [joint_pos_L(3,i) joint_pos_L(3,i)],'y')
%         plot3([joint_pos_R(1,i) joint_pos_R(1,i)], [joint_pos_R(2,i) joint_pos_R(2,i)+wing_l], [joint_pos_R(3,i) joint_pos_R(3,i)],'y')
%         plot3([-1.3671 -1.3671+BL], [0 0],[0 0],'y')
%         axis equal
%         hold off
%         
%         pause
        
        
    
    % Correction of the wing quaternions:
    for j = start:stop
       
        temp_qL = [pathDB.qL1(j,i) pathDB.qL2(j,i) pathDB.qL3(j,i) pathDB.qL4(j,i)];
        temp_qR = [pathDB.qR1(j,i) pathDB.qR2(j,i) pathDB.qR3(j,i) pathDB.qR4(j,i)];
        
        temp_qL = temp_qL./norm(temp_qL);
        temp_qR = temp_qR./norm(temp_qR);
        
        qL1_filt1(j,i) = temp_qL(1);
        qL2_filt1(j,i) = temp_qL(2);
        qL3_filt1(j,i) = temp_qL(3);
        qL4_filt1(j,i) = temp_qL(4);

        qR1_filt1(j,i) = temp_qR(1);
        qR2_filt1(j,i) = temp_qR(2);
        qR3_filt1(j,i) = temp_qR(3);
        qR4_filt1(j,i) = temp_qR(4);

        
        
        % Calculate the angle of attack and the sideslip angle w.r.t. the
        % body.
        
        U_vect = [pathDB.u_filt(j,i); pathDB.v_filt(j,i); pathDB.w_filt(j,i)];
        
        a_vect = [pathDB.ax_filt(j,i); pathDB.ay_filt(j,i); pathDB.az_filt(j,i)];
    
        DCM = quat2matNEW([pathDB.qb1_filt(j,i) pathDB.qb2_filt(j,i) pathDB.qb3_filt(j,i) pathDB.qb4_filt(j,i)]);
        
        U_body = DCM'*U_vect;
        
        a_body = DCM'*a_vect;
        
        alfa_body(j,i) = real(atan2(U_body(3),U_body(1)));
        
%         beta_body(j,i) = real(atan2(U_body(2),U_body(1)));
        
%         beta_body(j,i) = real(atan2(U_body(2),sqrt(U_body(1)^2+U_body(2)^2+U_body(3)^2)));

        beta_body(j,i) = real(asin(U_body(2)/sqrt(U_body(1)^2+U_body(2)^2+U_body(3)^2)));

        clear U_vect DCM
        
        u_body(j,i) = U_body(1);
        v_body(j,i) = U_body(2);
        w_body(j,i) = U_body(3);
        
        ax_body(j,i) = a_body(1);
        ay_body(j,i) = a_body(2);
        az_body(j,i) = a_body(3);
        
        
        % Calculate the roll, pitch and yaw angle of the body in the global
        % frame.

        
        roll_body(j,i) = atan2(2*(pathDB.qb4_filt(j,i).*pathDB.qb1_filt(j,i)+pathDB.qb2_filt(j,i).*pathDB.qb3_filt(j,i)),(1-2*(pathDB.qb1_filt(j,i).^2+pathDB.qb2_filt(j,i).^2)));
        pitch_body(j,i) = asin(2*(pathDB.qb4_filt(j,i).*pathDB.qb2_filt(j,i)-pathDB.qb3_filt(j,i).*pathDB.qb1_filt(j,i)));
        yaw_body(j,i) = atan2(2*(pathDB.qb4_filt(j,i).*pathDB.qb3_filt(j,i)+pathDB.qb1_filt(j,i).*pathDB.qb2_filt(j,i)),1-2*(pathDB.qb2_filt(j,i).^2+pathDB.qb3_filt(j,i).^2));      
        
    end
end


    temp.qL1_filt1 = qL1_filt1;
    temp.qL2_filt1 = qL2_filt1;
    temp.qL3_filt1 = qL3_filt1;
    temp.qL4_filt1 = qL4_filt1;

    temp.qR1_filt1 = qR1_filt1;
    temp.qR2_filt1 = qR2_filt1;
    temp.qR3_filt1 = qR3_filt1;
    temp.qR4_filt1 = qR4_filt1;
    
    temp.alfa_body = alfa_body;
    temp.beta_body = beta_body;
    
    temp.roll_body = roll_body;
    temp.pitch_body = pitch_body;
    temp.yaw_body = yaw_body;
    
    temp.u_body = u_body;
    temp.v_body = v_body;
    temp.w_body = w_body;
    
    temp.ax_body = ax_body;
    temp.ay_body = ay_body;
    temp.az_body = az_body;
    
    temp.body_length = body_length;
    temp.wing_length = wing_length;
    temp.joint_pos_L = joint_pos_L;
    temp.joint_pos_R = joint_pos_R;




% save(savefile,'qL1_filt1','qL2_filt1','qL3_filt1','qL4_filt1','qR1_filt1','qR2_filt1','qR3_filt1','qR4_filt1','alfa_body','beta_body','roll_body','pitch_body','yaw_body','u_body','v_body','w_body','ax_body','ay_body','az_body','body_length','wing_length','joint_pos_L','joint_pos_R')