



pattern.x_num = 96;     % There are 96 pixel around the display (12x8)

pattern.y_num = 12;      % 16 frames of Y makes one complete cycle

pattern.num_panels = 48;    % This is the number of unique Panel IDs required.

pattern.gs_val = 3;     % This pattern will use 8 intensity levels

pattern.row_compression = 1; % size of each frame is 4x96, because of row compression.

Pats = zeros(4, 96, pattern.x_num, pattern.y_num);  %initializes the array with zeros



PatsL = zeros(4, 96, pattern.x_num, pattern.y_num);

PatsR = zeros(4, 96, pattern.x_num, pattern.y_num);

Pats(:, :, 1, 1) = repmat((repmat([0 1 2 4 6 7 7 6 4 2 1 0], 4,8)), 1, 1);

PatsL(:, :, 1, 1) = Pats(:, :, 1, 1);

PatsR(:, :, 1, 1) = Pats(:, :, 1, 1);



for k=2:96 

    Pats(:,:,k,1) = ShiftMatrix(Pats(:,:,k-1,1),1,'r','y');

    PatsL(:,:,k,1) = ShiftMatrix(Pats(:,:,k-1,1),1,'r','y');

    PatsR(:,:,k,1) = ShiftMatrix(Pats(:,:,k-1,1),1,'r','y');

end



for j=2:12 

    PatsL(:,:,:,j) = ShiftMatrix(PatsL(:,:,:,j-1),1,'l','y');

    PatsR(:,:,:,j) = ShiftMatrix(PatsR(:,:,:,j-1),1,'r','y');

end



for k=2:96

    Pats(:,1:(k-1),k,:) = PatsL(:,1:(k-1),k,:);

    Pats(:,k:96,k,:) = PatsR(:,k:96,k,:);

end





pattern.Pats = Pats; 		% put data in structure

% A = 1:48;              	% define panel structure vector

% pattern.Panel_map = flipud(reshape(A, 4, 12));

pattern.Panel_map = [12 8 4 11 7 3 10 6 2  9 5 1; 24 20 16 23 19 15 22 18 14 21 17 13; 36 32 28 35 31 27 34 30 26 33 29 25; 48 44 40 47 43 39 46 42 38 45 41 37];
A = 1:120;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, 24, 5)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);



directory_name = 'C:\Users\shazam\Desktop\Patterns';
directory_name = cd;

str = [directory_name '\Pattern_cos_grating_expansion'];  % name must begin with ?Pattern_?

save(str, 'pattern');



