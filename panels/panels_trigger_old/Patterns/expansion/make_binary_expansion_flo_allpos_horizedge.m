
n_row = 5;
n_colom = 24;

pattern.x_num = 8*n_colom;     % expansion starts at all positions
pattern.y_num = 16;      % 12 frames of Y makes one complete cycle
pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary
pattern.row_compression = 1; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = zeros(n_row, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros

Pats(:, :, 1, 1) = repmat((repmat([0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 ], n_row,n_colom*8/pattern.y_num)), 1, 1);
Pats(1:3, :, 1, 1) = ShiftMatrix(Pats(1:3,:,1,1),8,'l','y');

PatsL = zeros(n_row, n_colom*8, 1, pattern.y_num);
PatsR = zeros(n_row, n_colom*8, 1, pattern.y_num);

PatsL(:, :, 1, 1) = Pats(:, :, 1, 1);
PatsR(:, :, 1, 1) = Pats(:, :, 1, 1);

for j=2:pattern.y_num

    PatsL(:,:,:,j) = ShiftMatrix(PatsL(:,:,:,j-1),1,'l','y');
    PatsR(:,:,:,j) = ShiftMatrix(PatsR(:,:,:,j-1),1,'r','y');

end
Pats(:,1:end/2,1,:) = PatsL(:,1:end/2,:,:);
Pats(:,end/2+1:end,1,:) = PatsR(:,end/2+1:end,:,:);

for k=2:pattern.x_num 

    Pats(:,:,k,:) = ShiftMatrix(Pats(:,:,k-1,:),1,'r','y');

end


% for k=2:pattern.x_num
% 
%     Pats(:,1:(k-1),k,:) = PatsL(:,1:(k-1),k,:);
%     Pats(:,k:end,k,:) = PatsR(:,k:end,k,:);
% 
% end

pattern.Pats = Pats; 		% put data in structure

A = 1:n_colom*n_row;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, n_colom, n_row)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);



directory_name = cd;
str = [directory_name '\Pattern_binary_expansion_flo_allpos_horizedge'];  % name must begin with ?Pattern_?

save(str, 'pattern');



