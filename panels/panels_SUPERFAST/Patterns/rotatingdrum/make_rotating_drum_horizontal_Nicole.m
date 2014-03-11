n_row = 5;
n_colom = 24;
SF = 16;
fr_max = SF;

pattern.x_num = 1;     % expansion starts at all vertical positions
pattern.y_num = fr_max;      % amount of frames (Y) for one complete expansion
pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary
pattern.row_compression = 0; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = zeros(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros
% Pats = zeros(L,M,N,O)
% L = number of panel rows
% M = number of panel columns
% N = frames in the 'x' direction
% O = frames in the 'y' direction

% % Single stripe
% Pats(:, :, 1, 1) = [ones(4,6) zeros(4,188)];
% % two stripes, 90 degs apart
% Pats(:, :, 1, 2) = [ones(4,6) zeros(4,6) ones(4,66) zeros(4,6) ones(4,12)]; 
% % 3 stripes
% Pats(:, :, 1, 3) = repmat([ones(4, 26) zeros(4,6)], 1, 3);

% 8 stripes
Pats(:, :, 1, 1) = repmat([ones(8, n_colom*8) zeros(8, n_colom*8)], n_row, 1);

for j = 2:fr_max
    Pats(:,:,1,j) = ShiftMatrix(Pats(:,:,j-1,1),1,'d','y');
    
end

% Pats(:,:,j,2) = ShiftMatrix(Pats(:,:,j-1,2),1,'r','y');
%     Pats(:,:,j,3) = ShiftMatrix(Pats(:,:,j-1,3),1,'r','y');
%     Pats(:,:,j,4) = ShiftMatrix(Pats(:,:,j-1,4),1,'r','y');

pattern.Pats = Pats;
A = [97:120; 73:96; 49:72; 25:48; 1:24]; % panel ID tags for Florian's rig
pattern.Panel_map = fliplr(A); % panel ID tag matrix flipped

pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);

directory_name = 'E:\Dropbox\WORK\panels\panels_SUPERFAST\Patterns\rotatingdrum';
str = [directory_name '\Pattern_8_wide_HORIZ_stripes_48_Pan']
save(str, 'pattern');
