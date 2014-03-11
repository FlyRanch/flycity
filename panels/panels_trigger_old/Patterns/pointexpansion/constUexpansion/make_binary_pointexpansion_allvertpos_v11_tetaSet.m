clear
clc

addpath('/home/florian/Dropbox/WORK/panels/panels/controller')
addpath('/home/florian/Dropbox/WORK/panels/panels/functions')
addpath('/home/florian/Dropbox/WORK/panels/panels/IO_tools')
addpath('/home/florian/Dropbox/WORK/panels/panels/Pattern_tools')
addpath('/home/florian/Dropbox/WORK/panels/panels/functions')

cd('C:\Users\florian\Dropbox\WORK\panels\panels\Patterns\pointexpansion\constUexpansion')
%% input
teta_crit = 65  % [deg] max response angle from Bender&Dickinson 2006
t_crit = 0.032 % [sec] time of teta_crit
fps = 72    % max fps without frame dropping

teta_max_set = teta_crit
t_max_set = 5588/7500/2; % sample time

% panel data
panel_max = 5; % amount of panels that are ON
poe_angle_hor = 0; % horizontal angle of point of expansion
% poe_angle_hor = 180; % horizontal angle of point of expansion

n_row = 5;
n_colom = 24;

%% pattern speed settings
panel_size = 32; % mm
pixperm = 8/panel_size * 1000;
pixperdeg = 8*n_colom/360;
pix_crit = teta_crit * pixperdeg;
fr_crit = t_crit*fps;

pixperfr = (pix_crit) / (fr_crit+1) % compensated for non-zero size at t=0

pix_max = teta_max_set * pixperdeg;
fr_max = t_max_set*fps;

% horizontal position of point of expansion
if poe_angle_hor == 0;
    poe_pix_hor = 8*9 + 4; % 0 point in cali
elseif poe_angle_hor == 180;
    poe_pix_hor = 8*3 -4 ; % right turn (180 deg from 0 point in cali)
end

% build pattern in center (no edge effects)
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


%% loop pattern
R=0; % expansion circle radius in pixels
dR = pixperfr/2; % dR = half the pixels per frame (r=d/2)

for j=2:pattern.y_num
    
%   [i j]
    Pats(:,:,:,j) = Pats(:,:,:,j-1);
    
    if 2*R <= pix_max
    R=R+dR; 
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

    % make time series db
    exp_steps(j-1,:) = [floor(2*R) floor(2*R)/pixperm floor(2*R)/pixperdeg (j-2)/fps];
end

%% shift to zero position
for k=1:pattern.x_num
    for l=1:pattern.y_num
        Pts(:,:) = Pats(:,:,k,l);
        Pts = ShiftMatrix(Pts,poe_pix_hor,'l','y');
        Pats(:,:,k,l) = Pts(:,:);
    end
end

%% make pattern structure
pattern.Pats = Pats; 		% put data in structure

A = 1:n_colom*n_row;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, n_colom, n_row)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);

%% settings dataset

fr_max = length(exp_steps);
t_max = exp_steps(end,4);

pix_max = min(8*panel_max,exp_steps(end,1));
teta_max = pix_max/pixperdeg;

settings.teta_crit = teta_crit;
settings.teta_max_set = teta_max_set;
settings.teta_max = teta_max;
settings.fps = fps;
settings.t_crit = t_crit;
settings.t_max = t_max;
settings.t_max_set = t_max_set;
settings.poe_angle_hor = poe_angle_hor;
settings.poe_pix_hor = poe_pix_hor;
settings.panel_max = panel_max;

%% save data
filename = ['pointexpansion_shiftTo' num2str(poe_angle_hor) 'deg_' num2str(round(teta_crit)) 'degAt' num2str(round(1000*t_crit)) 'ms_max' num2str(round(teta_max)) 'deg_tmax' num2str(round(1000*t_max)) 'ms_' num2str(fps) 'fps.mat'];
exp_variables = {'pixels' 'diameter' 'angle' 't'};

directory_name = cd;
str = [directory_name '/Pattern_' filename];  % name must begin with ?Pattern_?
save(str, 'pattern');

str = [directory_name '/expansion_timeseries_' filename];  %
save(str, 'settings', 'exp_variables', 'exp_steps');



