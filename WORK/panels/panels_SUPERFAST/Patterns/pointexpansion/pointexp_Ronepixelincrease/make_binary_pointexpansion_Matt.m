%% Add Paths
clear
clc

% Windows
% addpath('E:\Dropbox\WORK\panels\panels_SUPERFAST\controller')
% addpath('E:\Dropbox\WORK\panels\panels_SUPERFAST\functions')
% addpath('E:\Dropbox\WORK\panels\panels_SUPERFAST\IO_tools')
% addpath('E:\Dropbox\WORK\panels\panels_SUPERFAST\Pattern_tools')
% addpath('E:\Dropbox\WORK\panels\panels_SUPERFAST\functions')
% 
% cd('E:\Dropbox\WORK\panels\panels_SUPERFAST\Patterns\pointexpansion\pointexp_Ronepixelincrease')

% Ubuntu
addpath('/home/matt/Dropbox/WORK/panels/panels_SUPERFAST/Patterns/controller')
addpath('/home/matt/Dropbox/WORK/panels/panels_SUPERFAST/Patterns/functions')
addpath('/home/matt/Dropbox/WORK/panels/panels_SUPERFAST/Patterns/IO_tools')
addpath('/home/matt/Dropbox/WORK/panels/panels_SUPERFAST/Patterns/Pattern_tools')
addpath('/home/matt/Dropbox/WORK/panels/panels_SUPERFAST/Patterns/functions')

cd('/home/matt/Dropbox/WORK/panels/panels_SUPERFAST/Patterns/pointexpansion/pointexp_Ronepixelincrease')
%% input
fr_max = 99; % in BW max 100 frames fit on a panel

% panel data
n_row = 5;
n_colom = 24;

panel_size = 32; % mm
pixel = 8; % number of pixels per axis
pixperm = 8/panel_size * 1000;
pixperdeg = 8*n_colom/360;

mid_pos = 20; % vertical mid position

% horizontal angle of point of expansion
poe_angle_hor = 0; % 
% poe_angle_hor = 90; %
% poe_angle_hor = 180; %
% poe_angle_hor = 270; %

pos_head = 76; % pixel position that the fly is headed toward (pixel position of 0)

% horizontal position of point of expansion
if poe_angle_hor == 0;
    poe_pix_hor = pos_head; % 0 point in cali
elseif poe_angle_hor == 90;
    poe_pix_hor = (pos_head)-((n_colom* pixel)/4); % right turn (90 deg from 0 point in cali)
elseif poe_angle_hor == 180;
    poe_pix_hor = (pos_head)+(2*(n_colom * pixel)/4); % 180
elseif poe_angle_hor == 270;
    poe_pix_hor = (pos_head)+((n_colom * pixel)/4); % left turn (90 deg from 0 point in cali)
end

% build pattern in center (no edge effects)
poe_start = 8*n_colom/2; % middle of n coloms

pattern.x_num = 1;     % expansion starts at mid position
% pattern.x_num = 8*n_row;     % expansion starts at all vertical positions
pattern.y_num = fr_max;      % amount of frames (Y) for one complete expansion

pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary

% no row compression possible
% pattern.row_compression = 1; % size of each frame is 8*n_colom * 1*n_row, because of row compression.
pattern.row_compression = 0; % size of each frame is 8*n_colom * 8*n_row, because of NO row compression.

Pats = ones(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros


%% loop pattern
R=0; % expansion circle radius in pixels
dR = 1; % dR = ONE pixel per frame

for j=2:pattern.y_num
    j
    Pats(:,:,:,j) = Pats(:,:,:,j-1);
    R=R+dR; 
    
    for r=R-dR:R
        for teta_steps=0:10*360
            teta = teta_steps/10;
            
            a = round (mid_pos + r*cosd(teta));
            b = round (poe_start + r*sind(teta));
            [a b];
            
            if a>0 && a<=8*n_row && b>0 && b<=8*n_colom
                Pats(a,b,:,j) = 0;
            end
        end
    end

    % make time series db
    exp_steps(j-1,:) = [floor(2*R) floor(2*R)/pixperm floor(2*R)/pixperdeg];
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
settings.poe_angle_hor = poe_angle_hor;
settings.poe_pix_hor = poe_pix_hor;

%% save data
filename = ['pointexp_funcfr_pos' num2str(poe_angle_hor) 'deg_' num2str(round(fr_max)) 'fr_dR' num2str(dR) 'pix.mat'];
exp_variables = {'pixels' 'diameter' 'angle'};

directory_name = cd;
str = [directory_name '/Pattern_' filename];  % name must begin with ?Pattern_?
save(str, 'pattern');

str = [directory_name '/expansion_frames_' filename];  %
save(str, 'settings', 'exp_variables', 'exp_steps');



