% make_expanding_rectangle.m
% Make a pattern of a rectangle that will expand from the lower edge of the
% panel array. 

% Max Sizemore    September, 2011

pattern.x_num = 96; 		% (8*12) number of locations for the center of expansion
pattern.y_num = 32; 	% number of vertical positions 
pattern.num_panels = 48; 	% This is the number of unique Panel IDs required.
pattern.gs_val = 3; 	% This pattern will use 8 intensity levels
pattern.row_compression = 0;    % no row compression

StartWidth = 1; % starting width of the rectangle = 2*StartWidth + 2
StL = (48+1)-StartWidth; % Left start position
StR = 48+StartWidth; % Right start position

Pats = zeros(32, 96, pattern.x_num, pattern.y_num); % Initialize an array of zeros

Pats(32:32, (StL-1):(StR+1), 1, 1) = 1;
Pats(31:32, (StL-2):(StR+2), 1, 2) = 1;
Pats(30:32, (StL-3):(StR+3), 1, 3) = 1;
Pats(29:32, (StL-4):(StR+4), 1, 4) = 1;

Pats(32:32, (StL-1):(StR+1), 1, 2) = 3; %Pats(32:32, (StL-1):(StR+1), 1, 2) + 2;


Pats(31:32, (StL-2):(StR+2), 1, 3) = 3; %Pats(31:32, (StL-2):(StR+2), 1, 3) + 2;
Pats(32:32, (StL-1):(StR+1), 1, 3) = 5; %Pats(32:32, (StL-1):(StR+1), 1, 3) + 4;


Pats(30:32, (StL-3):(StR+3), 1, 4) = 3; %Pats(30:32, (StL-3):(StR+3), 1, 4) + 2;
Pats(31:32, (StL-2):(StR+2), 1, 4) = 5; %Pats(31:32, (StL-2):(StR+2), 1, 4) + 4;
Pats(32:32, (StL-1):(StR+1), 1, 4) = 7; %Pats(32:32, (StL-1):(StR+1), 1, 4) + 4;

for i = 5:32
        Pats((33-i):32, (StL-i):(StR+i), 1, i:32) = 1;
        Pats((34-i):32, (StL-i+1):(StR+i-1), 1, i) = 3;%Pats((34-i):32, (StL-i-5):(StR+i+5), 1, i) + 2;
        Pats((35-i):32, (StL-i+2):(StR+i-2), 1, i) = 5;%Pats((35-i):32, (StL-i-4):(StR+i+4), 1, i) + 2;
        Pats((36-i):32, (StL-i+3):(StR+i-3), 1, i) = 7;%Pats((36-i):32, (StL-i-3):(StR+i+3), 1, i) + 2;
%         Pats((37-i):32, (StL-i-2):(StR+i+2), 1, i) = 7;
end

for x = 2:96
    Pats(:, :, x, :) = ShiftMatrix(Pats(:, :, x-1, :), 1, 'r', 'y');
end

for x = 1:48
    Pats(:, 1:x, x, :)= 0;
end

for x = 49:96
    Pats(:, x:96, x, :) = 0;
end

pattern.Pats = Pats;

pattern.Panel_map = [12 8 4 11 7 3 10 6 2  9 5 1; 24 20 16 23 19 15 22 18 14 21 17 13; 36 32 28 35 31 27 34 30 26 33 29 25; 48 44 40 47 43 39 46 42 38 45 41 37];
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = make_pattern_vector(pattern);

directory_name = 'C:\Users\shazam\Desktop\Patterns';
str = [directory_name '\pattern_exp_rect_nowrap_graded'];
save(str, 'pattern');
