clear
clc

cd('C:\Users\florian\Dropbox\WORK\panels\panels_trigger\Patterns\test')

n_row = 5;
n_colom = 24;

pattern.x_num = 11;     % 1 frame on, 10 frames off
pattern.y_num = 1;      % no y variable
pattern.num_panels = n_colom*n_row;    % =24*5: This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will use 1 bit

Pats = zeros(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros
% make pattern:  1 frame off, 1 frame on, 10 frames off
Pats(:, :, 1, :) = ones(n_row*8, n_colom*8, 1, pattern.y_num);

pattern.Pats = Pats; 		% put data in structure

A = 1:120;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, 24, 5)',2);

pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);


directory_name = cd;
str = [directory_name '\Pattern_1on_10off'];  % name must begin with ?Pattern_?

save(str, 'pattern');



