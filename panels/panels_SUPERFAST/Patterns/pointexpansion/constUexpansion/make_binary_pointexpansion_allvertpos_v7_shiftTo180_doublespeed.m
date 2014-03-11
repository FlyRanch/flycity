clear
clc

% Bender&Dickinson 2006
teta_crit = 62  % [deg] max response angle
delta = .049 % [sec] delay: NOT USED!

% panel data
%% horizontal position of point of expansion 
% poe_hor = 8*9 + 4; % left turn (0 point in cali)
poe_hor = 8*3 -4 ; % right turn (180 deg from 0 point in cali)
% poe_hor = 8*n_colom/2; % panel mid
%%
fr_crit = 5; % teta_crit within 5 frames
fr_crit = 2.5; % teta_crit within 2.5 frames: DOUBLE SPEED!

fr_max = 4* fr_crit; 
d = 0.120; % [m] panel-fly distance
pixperm = 8/32 * 1000; % pixel size

R_crit = d * tand(teta_crit/2); % [m] critical radius
pix_crit = R_crit * pixperm;
pixperfr = pix_crit / fr_crit;

n_row = 5;
n_colom = 24;

% pattern.x_num = 1;     % expansion starts at mid position
pattern.x_num = 8*n_row;     % expansion starts at all vertical positions
pattern.y_num = fr_max;      % 20 frames of Y makes one complete expansion

pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary

% no row compression possible
% pattern.row_compression = 1; % size of each frame is 8*n_colom * 1*n_row, because of row compression.
pattern.row_compression = 0; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = ones(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros

%vertical mid point-of-expansion 
poe = 8*n_colom/2; % middle of 5 panels

for i=1:pattern.x_num
    R=-1; % expansion circle radius
    for j=2:pattern.y_num
        [i j];
        R=R+pixperfr;
        Pats(:,:,i,j) = Pats(:,:,i,j-1);
        for r=R-pixperfr:R
            for teta=0:1:360
                a = round (i + r*cosd(teta));
                b = round (poe + r*sind(teta));
                [a b];
                if a>0 && a<=8*n_row && b>0 && b<=8*n_colom
                    Pats(a,b,i,j) = 0;
                end
            end
        end
    end
end


%% shift to zero position

% % custom shift
% Pats_bkp = Pats;
% Pats(:,1:end-poe_hor,:,:) = Pats_bkp(:,poe_hor+1:end,:,:);
% Pats(:,end-poe_hor+1:end,:,:) = Pats_bkp(:,1:poe_hor,:,:);
    
for k=1:pattern.x_num
    for l=1:pattern.y_num
        Pts(:,:) = Pats(:,:,k,l);
%         Pts = ShiftMatrix(Pts,poe_hor,'l','y');
        Pts = ShiftMatrix(Pts,poe_hor,'r','y');     % 180deg from 0 cali
        Pats(:,:,k,l) = Pts(:,:);
    end
end

%% save data
pattern.Pats = Pats; 		% put data in structure

A = 1:n_colom*n_row;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, n_colom, n_row)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);



directory_name = cd;
str = [directory_name '\Pattern_binary_pointexpansion_allvertpos_shift180deg_DOUBLESPEED'];  % name must begin with ?Pattern_?

save(str, 'pattern');



