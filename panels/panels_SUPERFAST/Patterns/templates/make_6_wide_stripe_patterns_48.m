% make_6_wide_med_cont_pattern_48.m

pattern.x_num = 96; 	% There are 96 pixel around the display (12x8) 
pattern.y_num = 4; 		% two frames of Y, at 2 different spatial frequencies
pattern.num_panels = 48; 	% This is the number of unique Panel IDs required.
pattern.gs_val = 1; 	% This pattern will use 8 intensity levels
pattern.row_compression = 1;

Pats = zeros(4, 96, pattern.x_num, pattern.y_num);
Pats(:, :, 1, 1) = [ones(4,90) zeros(4,6)]; % one stripe
% two stripes, 90 degs apart
Pats(:, :, 1, 2) = [ones(4,6) zeros(4,6) ones(4,66) zeros(4,6) ones(4,12)]; 
% 3 stripes
Pats(:, :, 1, 3) = repmat([ones(4, 26) zeros(4,6)], 1, 3);
% 8 stripes
Pats(:, :, 1, 4) = repmat([ones(4, 6) zeros(4,6)], 1, 8);

for j = 2:96
    Pats(:,:,j,1) = ShiftMatrix(Pats(:,:,j-1,1),1,'r','y');
    Pats(:,:,j,2) = ShiftMatrix(Pats(:,:,j-1,2),1,'r','y');
    Pats(:,:,j,3) = ShiftMatrix(Pats(:,:,j-1,3),1,'r','y');
    Pats(:,:,j,4) = ShiftMatrix(Pats(:,:,j-1,4),1,'r','y');
end

pattern.Pats = Pats;
A = [1:12; 13:24; 25:36; 37:48]; % panel ID tags for Florian's rig
pattern.Panel_map = fliplr(A); % panel ID tag matrix flipped

% directory_name = 'C:\temp';
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);

directory_name = 'E:\Dropbox\WORK\panels\panels_SUPERFAST\Patterns\rotatingdrum';
str = [directory_name '\Pattern_6_wide_stripes_48_Pan']
save(str, 'pattern');
