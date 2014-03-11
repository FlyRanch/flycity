% makes a sine wave grating pattern with grayscale and row compression
% Max 9/2011

pattern.x_num = 96; 	% There are 96 pixel around the display (12x8) 
pattern.y_num = 6; 		% three frames of Y, each at 2 different bar widths
pattern.num_panels = 48; 	% This is the number of unique Panel IDs required.
pattern.gs_val = 3; 	% This pattern will use 8 intensity levels
pattern.row_compression = 1; % Row compression is on

nar = repmat([1 3 5 7 7 5 3 1],[4 1]);
wid = repmat([1 1 3 3 5 5 7 7 7 7 5 5 3 3 1 1],[4 1]);

Pats = zeros(4, 96, pattern.x_num, pattern.y_num); 	%initializes the array with zeros
% make grating patterns, periods are 120 and 60 degrees, using all 8 gscale values
Pats(:, 45:52, 1, 1) = nar;
Pats(:, 84:91, 1, 2) = nar;
Pats(:, 5:12, 1, 3) = nar;

Pats(:, 41:56, 1, 4) = wid;
Pats(:, 80:95, 1, 5) = wid;
Pats(:, 1:16, 1, 6) = wid;

for j = 2:96 			%use ShiftMatrixPats to rotate stripe image
    Pats(:,:,j,1) = ShiftMatrix(Pats(:,:,j-1,1),1,'r','y');
    Pats(:,:,j,2) = ShiftMatrix(Pats(:,:,j-1,2),1,'r','y');
    Pats(:,:,j,3) = ShiftMatrix(Pats(:,:,j-1,3),1,'r','y');
    Pats(:,:,j,4) = ShiftMatrix(Pats(:,:,j-1,4),1,'r','y');
    Pats(:,:,j,5) = ShiftMatrix(Pats(:,:,j-1,5),1,'r','y');
    Pats(:,:,j,6) = ShiftMatrix(Pats(:,:,j-1,6),1,'r','y');
end

Pats(:, 1:33, :, 2) = 0;
Pats(:, 53:96, :, 3) = 0;

Pats(:, 1:33, :, 5) = 0;
Pats(:, 53:96, :, 6) = 0;



pattern.Pats = Pats; 		% put data in structure 
% A = 1:48;              	% define panel structure vector
% pattern.Panel_map = flipud(reshape(A, 4, 12));
pattern.Panel_map = [12 8 4 11 7 3 10 6 2  9 5 1; 24 20 16 23 19 15 22 18 14 21 17 13; 36 32 28 35 31 27 34 30 26 33 29 25; 48 44 40 47 43 39 46 42 38 45 41 37];
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = make_pattern_vector(pattern);
directory_name = 'C:\Users\shazam\Desktop\Patterns';
str = [directory_name '\Pattern_vertbar_graded']; 	% name must begin with ‘Pattern_’
save(str, 'pattern');
