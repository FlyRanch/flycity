% function body_and_wing_model_v2_NOkin( settings )
clear
clc

    load('flyVars_addedmassDisc.mat')
    savefile = 'bodyNwingModel_4qsModel.mat';

    % Rotation matrices of body and wing reference frame and strokeplane
    % reference frame will be given. Body and wing mass and
    % inertia will be computed.
    

    % constants in correct quantities
    rho_air         = rho_air*1e-9;
    rho_cuticle     = rho_cuticle*1e-9;
    g               = g*1000;
    beta_strk       = -(body_angle/180)*pi;
    w_length = Lwing*1000;  % [mm] mean wing-length all wingbeats

    % Compute Rotation Matrix:
    Rstr = [ cos(beta_strk) 0 -sin(beta_strk); ...
                     0              1              0; ...
                     sin(beta_strk) 0 cos(beta_strk)];
    
    % Compute dimensions, mass, cg and inertia of body and wing:
    
    body_model.mass_fly           = nan(N,1);
    body_model.mass_body          = nan(N,1);
    body_model.length             = nan(N,1);
    body_model.Inertia            = nan(3,3,N);
    body_model.Joint_left         = nan(N,3);
    body_model.Joint_right        = nan(N,3);
    body_model.cg                 = nan(N,3);
    body_model.x_mod              = nan(21,13,N);
    body_model.y_mod              = nan(21,13,N);
    body_model.z_mod              = nan(21,13,N);
    
    wing_model.length             = nan(N,1);
    wing_model.mass               = nan(N,1);
    wing_model.virtual_mass       = nan(N,1);
    wing_model.Inertia            = nan(3,3,N);
    wing_model.virtual_Inertia    = nan(3,3,N);
    wing_model.area               = nan(N,1);
    
    wing_model.wing_cg_L          = nan(N,3);
    wing_model.y_sect_L           = nan(nr_chord_sect,3,N);
    wing_model.chords_L           = nan(N,nr_chord_sect);
    wing_model.x_LE_L             = nan(N,nr_chord_sect);
    wing_model.x_mod_L            = nan(25,2,N);
    wing_model.y_mod_L            = nan(25,2,N);
    wing_model.z_mod_L            = nan(25,2,N);
    
    wing_model.wing_cg_R          = nan(N,3);
    wing_model.y_sect_R           = nan(nr_chord_sect,3,N);
    wing_model.chords_R           = nan(N,nr_chord_sect);
    wing_model.x_LE_R             = nan(N,nr_chord_sect);
    wing_model.x_mod_R            = nan(25,2,N);
    wing_model.y_mod_R            = nan(25,2,N);
    wing_model.z_mod_R            = nan(25,2,N);
    
        %% Load body model:        
        
        xh = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1];
        
        [x_mod,y_mod,z_mod,mod_fit] = load_body_model_basic2( xh );
        
%% scale to average winglength
y_mod_R = y_mod{3};
w_length_basic           = max(y_mod_R(:,1))-min(y_mod_R(:,1));
Rwing = w_length/w_length_basic;

body_model.x_mod = x_mod{1} * Rwing;    % [mm]
body_model.y_mod = y_mod{1} * Rwing;
body_model.z_mod = z_mod{1} * Rwing;

wing_model.x_mod_L = x_mod{2} * Rwing;
wing_model.y_mod_L = y_mod{2} * Rwing;
wing_model.z_mod_L = z_mod{2} * Rwing;

wing_model.x_mod_R = x_mod{3} * Rwing;
wing_model.y_mod_R = y_mod{3} * Rwing;
wing_model.z_mod_R = z_mod{3} * Rwing;

body_model.length           = max(body_model.x_mod(:,1))-min(body_model.x_mod(:,1));
wing_model.length           = max(wing_model.y_mod_R(:,1))-min(wing_model.y_mod_R(:,1));

body_model.Joint_left       = body_model.length.*([0.2021 -0.1055 -0.1477])';
body_model.Joint_right      = body_model.length.*([0.2021 0.1055 -0.1477])';

djoint = 5.5/23 * wing_model.length;
body_model.Joint_left(2)        = -djoint;
body_model.Joint_right(2)       = djoint;

% body_model.mass_fly        = C_mfly*(wing_model.length/3)^3; % [kg]
body_model.mass_fly        = Mfly; % [kg]


%% inertia
        body_model_temp.x_mod = body_model.x_mod;
        body_model_temp.y_mod = body_model.y_mod;
        body_model_temp.z_mod = body_model.z_mod;
        wing_model_temp.x_mod = wing_model.x_mod_R;
        wing_model_temp.y_mod = wing_model.y_mod_R;
        wing_model_temp.z_mod = wing_model.z_mod_R;
        
%         figure()
%         hold on
%         surf(body_model_temp.x_mod,body_model_temp.y_mod,body_model_temp.z_mod,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
% %         surf(body_model.Joint_left(i,1)+wing_model_temp.x_mod,body_model.Joint_left(i,2)+wing_model_temp.y_mod,body_model.Joint_left(i,3)+wing_model_temp.z_mod,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
%         surf(wing_model_temp.x_mod,wing_model_temp.y_mod,wing_model_temp.z_mod,'facecolor',[0.3 0.3 0.3],'edgecolor','k','facelighting','phong');
%         plot3(body_model.Joint_left(i,1),body_model.Joint_left(i,2),body_model.Joint_left(i,3),'o')
%         axis equal
%         hold off
%         
%         pause
        
        [cg_body,cg_L,cg_R,I_body,I_wing,I_v_wing,m_w,m_v_w,y_sect,chords,area_w,x_LE] = ...
            comp_Inertia(rho_air,rho_cuticle,h_wing,body_model.length,wing_model.length,body_model.mass_fly, ...
            body_model_temp,wing_model_temp,body_model.Joint_right,nr_chord_sect);
        
%         figure()
%         hold on
%         surf(x_mod{1},y_mod{1},z_mod{1},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong')
%         surf(x_mod{2},y_mod{2},z_mod{2},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong')
%         surf(x_mod{3},y_mod{3},z_mod{3},'facecolor',[0.5 0.5 0.5],'edgecolor','k','facelighting','phong')
%         axis equal
%         hold off

        body_model.cg                   = cg_body;
        body_model.Inertia            = I_body;
        body_model.mass_body              = body_model.mass_fly-2*m_w;

        wing_model.mass                   = m_w;
        wing_model.virtual_mass           = m_v_w;
        wing_model.Inertia            = I_wing;
        wing_model.virtual_Inertia    = I_v_wing;
        wing_model.area                   = area_w;
        wing_model.wing_cg_L            = cg_L;
        wing_model.wing_cg_R            = cg_R;
        wing_model.y_sect_L           = -y_sect';
        wing_model.chords_L             = chords;
        wing_model.x_LE_L               = x_LE;
        wing_model.y_sect_R           = y_sect';
        wing_model.chords_R             = chords;
        wing_model.x_LE_R               = x_LE;
        
        settings.rho_air = rho_air;
%         settings.C_mfly = C_mfly;
        settings.rho_cuticle = rho_cuticle;
        settings.h_wing = h_wing;
        settings.nr_chord_sect = nr_chord_sect;
        settings.g = g;
        settings.beta_strk = beta_strk;
        settings.Rstr = Rstr;
        
    % Save the structures:
    save(savefile,'body_model','wing_model','settings')


