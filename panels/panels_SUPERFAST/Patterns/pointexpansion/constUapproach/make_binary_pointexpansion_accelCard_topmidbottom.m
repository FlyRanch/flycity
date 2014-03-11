clear
clc

pixpermm = 8/32;

% virtual object (from Card & Dickinson 2008)
R = 35; %mm
D0 = 420; %mm
Dend = 420-250; %mm
dt = 300; %ms
U = (D0-Dend)/dt; %m/s

t=[0:dt/20:2*dt]'; %ms
% time until contact
t=[0:dt/20:D0/U]'; %ms

% panel projection
d = 120; %mm
h0 = 2.5*32; %mm

r = d*R./(D0 - U*t); %mm
h = h0*(D0-U*t)/D0; %mm

rpix = r .* pixpermm;
hpix = h .* pixpermm;

% panels
n_row = 5;
n_colom = 24;

% pattern.x_num = 1;     % expansion starts at mid position
% pattern.x_num = 8*n_row;     % expansion starts at all vertical positions
% pattern.y_num = 20;      % 20 frames of Y makes one complete expansion

% pattern.x_num = 1;     % expansion starts at middle position
pattern.x_num = 3;     % expansion starts at bottom, middle & top position
pattern.y_num = length(t)+1;      % t+1 frames of Y makes one complete expansion (no pattern at 1st frame)

pattern.num_panels = n_colom*n_row;    % This is the number of unique Panel IDs required.
pattern.gs_val = 1;     % This pattern will be binary

% no row compression possible
% pattern.row_compression = 1; % size of each frame is 8*n_colom * 1*n_row, because of row compression.
pattern.row_compression = 0; % size of each frame is 8*n_colom * 1*n_row, because of row compression.

Pats = ones(n_row*8, n_colom*8, pattern.x_num, pattern.y_num);  %initializes the array with zeros

%horizontal mid point-of-expansion 
poe_hor = 8*n_colom/2; % panel mid
%vertical mid point-of-expansion 
poe_ver = 8*n_row/2; % panel mid

r_prev = 0; % previous radius (for first step)

%% top expansion
for i=1
    for j=2:pattern.y_num
        [i j];
%         Pats(:,:,i,j) = Pats(:,:,i,j-1);
        for rnow=0:rpix(j-1)
            for teta=0:1:360
                a = round (poe_ver - hpix(j-1) + rnow*cosd(teta));
                b = round (poe_hor + rnow*sind(teta));
                [a b]
                if a>0 && a<8*n_row && b>0 && b<8*n_colom
                    Pats(a,b,i,j) = 0;
                end
            end
        end
        r_prev = rpix(j-1);
    end
end

%% mid pos
for i=2
    for j=2:pattern.y_num
        [i j];
        Pats(:,:,i,j) = Pats(:,:,i,j-1);
        for rnow=0:rpix(j-1)
            for teta=0:1:360
                a = round (poe_ver + rnow*cosd(teta));
                b = round (poe_hor + rnow*sind(teta));
                [a b]
                if a>0 && a<8*n_row && b>0 && b<8*n_colom
                    Pats(a,b,i,j) = 0;
                end
            end
        end
        r_prev = rpix(j-1);
    end
end

%% bottom
for i=3
    for j=2:pattern.y_num
        [i j];
%         Pats(:,:,i,j) = Pats(:,:,i,j-1);
        for rnow=0:rpix(j-1)
            for teta=0:1:360
                a = round (poe_ver + hpix(j-1) + rnow*cosd(teta));
                b = round (poe_hor + rnow*sind(teta));
                [a b]
                if a>0 && a<8*n_row && b>0 && b<8*n_colom
                    Pats(a,b,i,j) = 0;
                end
            end
        end
        r_prev = rpix(j-1);
    end
end

%last frame black
Pats(:,:,:,j) = 0;

%% shift to zero position
% for k=size(Pats,2):1
% 
% Pats(:,k,1,:) = ShiftMatrix(Pats(:,k,1,:),poe_hor,'l',1);
% Pats(:,:,2,:) = ShiftMatrix(Pats(:,:,2,:),poe_hor,'l',1);
% Pats(:,:,3,:) = ShiftMatrix(Pats(:,:,3,:),poe_hor,'l',1);
% 
% end

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

pattern.Pats = pt; 		% put data in structure

A = 1:n_colom*n_row;              	% define panel structure vector
pattern.Panel_map = rot90(reshape(A, n_colom, n_row)',2);

pattern.BitMapIndex = process_panel_map(pattern);

pattern.data = Make_pattern_vector(pattern);



directory_name = cd;
str = [directory_name '\Pattern_binary_pointexpansion_accelCard_topmidbottom_shift'];  % name must begin with ?Pattern_?

save(str, 'pattern');



