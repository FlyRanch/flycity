clear
clc

n_row = 5;
n_colom = 24;

% pattern.x_num = 1;     % expansion starts at mid position
pattern.x_num = 8*n_row;     % expansion starts at all vertical positions
pattern.y_num = 20;      % 20 frames of Y makes one complete expansion

pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary

% no row compression possible
% pattern.row_compression = 1; % size of each frame is 8*n_colom * 1*n_row, because of row compression.
pattern.row_compression = 0; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = ones(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros

r=-1; % expansion circle radius
poe = [8*n_row/2 8*n_colom/2] % mid point-of-expansion

for j=2:pattern.y_num
    r=r+1;
    Pats(:,:,round(poe(1)),j) = Pats(:,:,round(poe(1)),j-1);
    
    for teta=0:360
        a = round (poe(1) + r*cosd(teta));
        b = round (poe(2) + r*sind(teta));
        if a>0 && a<8*n_row && b>0 && b<8*n_colom
            Pats(a,b,round(poe(1)),j) = 0;
        end
    end
end

% shift up
for k=round(poe(1))+1:pattern.x_num 

    Pats(:,:,k,:) = ShiftMatrix(Pats(:,:,k-1,:),1,'u',1);

end

% shift down
for k=round(poe(1))-1:-1:1

    Pats(:,:,k,:) = ShiftMatrix(Pats(:,:,k+1,:),1,'d',1);

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
str = [directory_name '\Pattern_binary_pointexpansion_allvertpos'];  % name must begin with ?Pattern_?

save(str, 'pattern');



