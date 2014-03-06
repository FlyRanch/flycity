



pattern.x_num = 192;     % =24*8: There are 24 panels x 8 leds around the display

pattern.y_num = 12;      % 12 frames of Y makes one complete cycle

pattern.num_panels = 24*5;    % This is the number of unique Panel IDs required.

pattern.gs_val = 3;     % This pattern will use 8 intensity levels

pattern.row_compression = 1; % size of each frame is 5x96, because of row compression.

Pats = zeros(5, 192, pattern.x_num, pattern.y_num);  %initializes the array with zeros



PatsL = zeros(5, 192, pattern.x_num, pattern.y_num);

PatsR = zeros(5, 192, pattern.x_num, pattern.y_num);

Pats(:, :, 1, 1) = repmat((repmat([0 1 2 4 6 7 7 6 4 2 1 0], 5,16)), 1, 1);

PatsL(:, :, 1, 1) = Pats(:, :, 1, 1);

PatsR(:, :, 1, 1) = Pats(:, :, 1, 1);



for k=2:size(PatsL,3) 

    Pats(:,:,k,1) = ShiftMatrix(Pats(:,:,k-1,1),1,'r','y');

    PatsL(:,:,k,1) = ShiftMatrix(Pats(:,:,k-1,1),1,'r','y');

    PatsR(:,:,k,1) = ShiftMatrix(Pats(:,:,k-1,1),1,'r','y');

end



for j=2:size(PatsL,4) 

    PatsL(:,:,:,j) = ShiftMatrix(PatsL(:,:,:,j-1),1,'l','y');

    PatsR(:,:,:,j) = ShiftMatrix(PatsR(:,:,:,j-1),1,'r','y');

end



for k=2:size(PatsL,3)

    Pats(:,1:(k-1),k,:) = PatsL(:,1:(k-1),k,:);

    Pats(:,k:end,k,:) = PatsR(:,k:end,k,:);

end





pattern.Pats = Pats; 		% put data in structure

A = 1:120;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, 24, 5)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);


directory_name = cd;
str = [directory_name '\Pattern_cos_grating_expansion'];  % name must begin with ?Pattern_?

save(str, 'pattern');



