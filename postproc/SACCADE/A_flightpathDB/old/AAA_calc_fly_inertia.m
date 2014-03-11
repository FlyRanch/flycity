% calc body & wing & added mass inertia
clear
clc

%% const
rho_water=1000;
rho_air=1.225;

% fly data
Lwing = .00299; %m
ARwing_fly  = 3.1579; % wing aspect ratio of fly (based on Elzinga's R/c)
Lbody = 3.138e-3;
body_angle = 47.5; 

% strokeplane
strokeplane_WBkin = -47.5;

% rotation axis


% 8 fem & 7 male, starved for 24h, by FTM. 50%male/female
g = 9.81
Mfly = 1.8e-6;
Mg_fly = Mfly*g;

rbody = sqrt(Mfly/rho_water/pi/Lbody);
cwing = Lwing/ARwing_fly;

% wing mass 1% of body mass
Mwing=0.005*Mfly;
Maddedmass = .25*pi*cwing^2*Lwing*rho_air;
MwingNair = Mwing+Maddedmass;

%% damping coeff [Nm/s] !!! RAD based !!!

% Yaw damping from Dickinson et al 2010 & converted to Dh scale
Cyaw   = 1.14e-10;

% roll and pitch damping scaled from yaw acordoing to Cheng et al 2009 (table 2)
Croll  = 0.55 * Cyaw;
Cpitch = 0.19 * Cyaw;

%% body inertia: model as cylinder at 47.5 degrees
% inertial moment coeff [Nm/s^2] !!! RAD based !!!

Ix_cyl = 1/12 * Mfly * (3*rbody^2 + Lbody^2);
Iz_cyl = 1/2 * Mfly * rbody^2;

Iyaw_body  = (Ix_cyl + Iz_cyl)/2 - ((Ix_cyl - Iz_cyl)/2)*cosd(180-2*body_angle)
Iroll_body = (Ix_cyl + Iz_cyl)/2 - ((Ix_cyl - Iz_cyl)/2)*cosd(2*body_angle)
Ipitch_body = Ix_cyl

%% wing inertia: modeled as wing disk

% single disk ring
Iyaw_wing = .5*MwingNair*( (Lwing+rbody)^2 + rbody^2);
Iroll_wing = .25*MwingNair*( (Lwing+rbody)^2 + rbody^2);
Ipitch_wing = Iroll_wing;

% N disk rings
N=100;
for n=1:N
    r_out = rbody + n/N * Lwing;
    r_in = rbody + (n-1)/N * Lwing;
    
    Iyaw_wing_sub(n,1)  =   .5*MwingNair/N*(r_out^2 + r_in^2);
    Iroll_wing_sub(n,1) =  .25*MwingNair/N*(r_out^2 + r_in^2);
    Ipitch_wing_sub(n,1) = .25*MwingNair/N*(r_out^2 + r_in^2);
end

Iyaw_wing = sum(Iyaw_wing_sub)
Iroll_wing = sum(Iroll_wing_sub)
Ipitch_wing = sum(Ipitch_wing_sub)

%% added mass modeled as a sphere minus cylinder
Maddedsphere = rho_air * 4/3 * pi * (Lwing + rbody)^3;
Iaddedsphere = 2/5* Maddedsphere * (Lwing + rbody)^2;

Mairbody = pi*rbody^2*Lbody*rho_air;
Ix_airbody = 1/12 * Mairbody * (3*rbody^2 + Lbody^2);
Iz_airbody = 1/2 * Mairbody * rbody^2;

Iyaw_airbody  = (Ix_airbody + Iz_airbody)/2 - ((Ix_airbody - Iz_airbody)/2)*cosd(180-2*body_angle);
Iroll_airbody = (Ix_airbody + Iz_airbody)/2 - ((Ix_airbody - Iz_airbody)/2)*cosd(2*body_angle);
Ipitch_airbody = Ix_airbody;


%% total inertia

% body & wing_addedmass disk
Iyaw_disc = Iyaw_body + Iyaw_wing
Iroll_disc = Iroll_body + Iroll_wing
Ipitch_disc = Ipitch_body + Ipitch_wing

% body & addedmass sphere-cylinder
Iyaw_sphere = Iyaw_body + Iaddedsphere - Iyaw_airbody
Iroll_sphere = Iroll_body + Iaddedsphere - Iroll_airbody
Ipitch_sphere = Ipitch_body + Iaddedsphere - Ipitch_airbody


%% save data
% body & wing_addedmass disk
Iyaw = Iyaw_disc;
Iroll = Iroll_disc;
Ipitch = Ipitch_disc;
save('flyVars_addedmassDisc.mat')

% sphere
Iyaw = Iyaw_sphere;
Iroll = Iroll_sphere;
Ipitch = Ipitch_sphere;
save('flyVars_addedmassSphere.mat')

% MEAN disc&sphere
Iyaw = mean([Iyaw_disc;Iyaw_sphere])
Iroll = mean([Iroll_disc;Iroll_sphere])
Ipitch = mean([Ipitch_disc;Ipitch_sphere])
save('flyVars_addedmassDiscSphere.mat')








