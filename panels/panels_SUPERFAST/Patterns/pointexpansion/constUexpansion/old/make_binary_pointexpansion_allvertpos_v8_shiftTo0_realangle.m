clear
clc

addpath('/home/florian/Dropbox/WORK/panels/panels/controller')
addpath('/home/florian/Dropbox/WORK/panels/panels/functions')
addpath('/home/florian/Dropbox/WORK/panels/panels/IO_tools')
addpath('/home/florian/Dropbox/WORK/panels/panels/Pattern_tools')
addpath('/home/florian/Dropbox/WORK/panels/panels/functions')


%% input
teta_crit = 60  % [deg] max response angle from Bender&Dickinson 2006
fps = 48    % max fps without frame dropping
t_crit = 0.1 % [sec] time of teta_crit
panel_max = 11; % amount of panels on
poe_angle_hor = 0; % horizontal angle of point of expansion
% poe_angle_hor = 180; % horizontal angle of point of expansion

%% pattern speed settings
fr_crit = round(fps*t_crit);
t_crit_real = fr_crit/fps
fr_max = 4* fr_crit; 
t_max = fr_max/fps;

%% panel data
n_row = 5;
n_colom = 24;

pix_max = 8*11; % max pattern size

% horizontal position of point of expansion
if poe_angle_hor == 0;
    poe_hor = 8*9 + 4; % 0 point in cali
elseif poe_angle_hor == 180;
    poe_hor = 8*3 -4 ; % right turn (180 deg from 0 point in cali)
end

panel_size = 32;
pixperm = 8/panel_size * 1000;
pixperdeg = 8*n_colom/360;
pix_crit = teta_crit * pixperdeg;
pixperfr = pix_crit / fr_crit;

% make pattern in center (no edge effects)
poe_start = 8*n_colom/2; % middle of n coloms

% pattern.x_num = 1;     % expansion starts at mid position
pattern.x_num = 8*n_row;     % expansion starts at all vertical positions
pattern.y_num = fr_max;      % amount of frames (Y) for one complete expansion

pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary

% no row compression possible
% pattern.row_compression = 1; % size of each frame is 8*n_colom * 1*n_row, because of row compression.
pattern.row_compression = 0; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = ones(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros

%% settings dataset
settings.teta_crit = teta_crit;
settings.fps = fps;
settings.t_crit = t_crit_real;
settings.t_max = t_max;
settings.poe_angle_hor = poe_angle_hor;
settings.poe_pix_hor = poe_hor;
settings.panel_max = panel_max;

exp_variables = {'pixels' 'diameter' 'angle' 't'};

%% loop pattern
dR = pixperfr/2;
R=0; % expansion circle radius in pixels
for j=2:pattern.y_num
%         [i j]
    R=R+dR; % update R with half the pixels per frame (r=d/2)

    % make time series db
    if 2*R<pix_max
        exp_steps(j-1,:) = [2*R 2*R/pixperm 2*R/pixperdeg (j-2)/fps];
    else
        exp_steps(j-1,:) = [pix_max pix_max/pixperm pix_max/pixperdeg (j-2)/fps];
    end

    Pats(:,:,:,j) = Pats(:,:,:,j-1);
    for r=R-dR:R
        for teta=0:1:360
            for i=1:pattern.x_num
                a = round (i + r*cosd(teta));
                b = round (poe_start + r*sind(teta));
                [a b];
                if a>0 && a<=8*n_row && b>0 && b<=8*n_colom
                    Pats(a,b,i,j) = 0;
                end
            end
        end
    end
end

%% shift to zero position
for k=1:pattern.x_num
    for l=1:pattern.y_num
        Pts(:,:) = Pats(:,:,k,l);
        Pts = ShiftMatrix(Pts,poe_hor,'l','y');
        Pats(:,:,k,l) = Pts(:,:);
    end
end

%% make pattern structure
pattern.Pats = Pats; 		% put data in structure

A = 1:n_colom*n_row;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, n_colom, n_row)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);

%% save data
directory_name = cd;
str = [directory_name '/Pattern_binary_pointexpansion_allvertpos_shiftTo0_realangle'];  % name must begin with ?Pattern_?
save(str, 'pattern');

str = [directory_name '/expansion_steps_binary_pointexpansion_allvertpos_shiftTo0_realangle'];  %
save(str, 'settings', 'exp_variables', 'exp_steps');


