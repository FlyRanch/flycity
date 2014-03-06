clear
clc

addpath('/home/florian/Dropbox/WORK/panels_SUPERFAST/panels/controller')
addpath('/home/florian/Dropbox/WORK/panels_SUPERFAST/panels/functions')
addpath('/home/florian/Dropbox/WORK/panels_SUPERFAST/panels/IO_tools')
addpath('/home/florian/Dropbox/WORK/panels_SUPERFAST/panels/Pattern_tools')
addpath('/home/florian/Dropbox/WORK/panels_SUPERFAST/panels/functions')

cd('C:\Users\florian\Dropbox\WORK\panels\panels_SUPERFAST\Patterns\test_setup')
%% input
n_row = 1;
n_colom = 1;

pattern.x_num = 1;     % single expansion
pattern.y_num = 9;      % amount of frames (Y) for one complete expansion

pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 4;     % This pattern will NOT be binary

% no row compression possible
pattern.row_compression = 0; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = 8*ones(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros


%% loop pattern
for j=2:pattern.y_num
    
%     Pats(:,:,:,j) = Pats(:,:,:,j-1);
    Pats(1:(j-1),1:(j-1),1,j) = 0;
end

%% make pattern structure
pattern.Pats = Pats; 		% put data in structure

A = 1:n_colom*n_row;              	% define panel structure vector
% A = 127;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, n_colom, n_row)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);

%% save data
filename = ['test_mirrorpanels_SUPERFAST'];

directory_name = cd;
str = [directory_name '/Pattern_' filename];  % name must begin with ?Pattern_?
save(str, 'pattern');
