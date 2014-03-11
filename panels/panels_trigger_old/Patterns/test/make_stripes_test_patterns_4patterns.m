% make_stripes_test_patterns.m

% 12 panel pattern
pattern.x_num = 4;
pattern.y_num = 1;
pattern.num_panels = 120;
pattern.gs_val = 1;

InitPat = repmat([zeros(5*8,2), ones(5*8,2)], 1,48);
Pats = zeros(5*8, 192, pattern.x_num, pattern.y_num);

Pats(:,:,1,1) = InitPat;

for j = 2:pattern.x_num
    Pats(:,:,j,1) = ShiftMatrix(Pats(:,:,j-1,1), 1, 'r', 'y'); 
end

pattern.Pats = Pats;

% pattern.Panel_map = [1 2 3 4 5 6 7 8 9 10 11 12];
A = 1:120;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, 24, 5)',2);

pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);

directory_name = cd;
str = [directory_name '\Pattern_stripes_4patterns']
save(str, 'pattern');

