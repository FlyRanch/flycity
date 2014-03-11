% makes a sine wave grating pattern with grayscale and row compression
% Max 9/2011
pattern.x_num = 96; 	% There are 96 pixel around the display (12x8) 
pattern.y_num = 2; 		% two frames of Y, at 2 different spatial frequencies
pattern.num_panels = 48; 	% This is the number of unique Panel IDs required.
pattern.gs_val = 3; 	% This pattern will use 8 intensity levels
pattern.row_compression = 1;

Pats = zeros(4, 96, pattern.x_num, pattern.y_num); 	%initializes the array with zeros
% make grating patterns, periods are 120 and 60 degrees, using all 8 gscale values
Pats(:, :, 1, 1) = repmat(round(3.5*(sin((4*pi/96)*[0:95])+1) ), 4, 1);
Pats(:, :, 1, 2) = repmat(round(3.5*(sin((8*pi/96)*[0:95])+1) ), 4, 1);

for j = 2:96 			%use ShiftMatrixPats to rotate stripe image
    Pats(:,:,j,1) = ShiftMatrix(Pats(:,:,j-1,1),1,'r','y');
    Pats(:,:,j,2) = ShiftMatrix(Pats(:,:,j-1,2),1,'r','y');
end

pattern.Pats = Pats; 		% put data in structure 
% A = 1:48;              	% define panel structure vector
% pattern.Panel_map = flipud(reshape(A, 4, 12));
pattern.Panel_map = [12 8 4 11 7 3 10 6 2  9 5 1; 24 20 16 23 19 15 22 18 14 21 17 13; 36 32 28 35 31 27 34 30 26 33 29 25; 48 44 40 47 43 39 46 42 38 45 41 37];
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = make_pattern_vector(pattern);
directory_name = 'C:\Users\shazam\Desktop\Patterns';
str = [directory_name '\Pattern_singrating_horiz'] 	% name must begin with ‘Pattern_’
save(str, 'pattern');
