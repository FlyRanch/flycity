clear
clc

pixpermm = 8/32;

% virtual object (from Card & Dickinson 2008)
R = 35; %mm
D0 = 420; %mm
Dend = 420-250; %mm
dt = 300; %ms
U = (D0-Dend)/dt; %m/s


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

Rstep=2; % expansion stepsize (amount of leds)
%horizontal mid point-of-expansion 
poe = 8*n_colom/2; % panel mid

for i=1:pattern.x_num
    R=-1; % expansion circle radius
    for j=2:pattern.y_num
        [i j];
        R=R+Rstep;
        Pats(:,:,i,j) = Pats(:,:,i,j-1);

        for r=R-Rstep:R
            for teta=0:1:360
                a = round (i + r*cosd(teta));
                b = round (poe + r*sind(teta));
                [a b]
                if a>0 && a<8*n_row && b>0 && b<8*n_colom
                    Pats(a,b,i,j) = 0;
                end
            end
        end
    end
end

% % shift up
% for k=round(poe(1))+1:pattern.x_num 
% 
%     Pats(:,:,k,:) = ShiftMatrix(Pats(:,:,k-1,:),1,'u',1);
% 
% end
% 
% % shift down
% for k=round(poe(1))-1:-1:1
% 
%     Pats(:,:,k,:) = ShiftMatrix(Pats(:,:,k+1,:),1,'d',1);
% 
% end


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



