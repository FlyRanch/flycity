% collision avoidance 
clear
clc

load('kine_20130122_S0006_bodywidth.mat')

b1 = (data.kine.FLYBODY1.data.length(data.kine.FLYBODY1.data.length~=0));
b2 = (data.kine.FLYBODY2.data.length(data.kine.FLYBODY2.data.length~=0));
b = mean([b1 b2]);


load('kine_20130122_S0006_2flies.mat')
t = [0:5587]'/7500;

% frames for length calc
frames1 = find(data.kine.FLYBODY1.data.length~=0)';
frames2 = find(data.kine.FLYBODY2.data.length~=0)';
l1 = mean(data.kine.FLYBODY1.data.length(frames1));
l2 = mean(data.kine.FLYBODY2.data.length(frames2));

% frames for pos calc
x1(:,:) = data.kine.FLYBODY1.data.coords(1,1,:);
x2(:,:) = data.kine.FLYBODY2.data.coords(1,1,:);
frames1 = find(x1~=0);
frames2 = find(x2~=0);

if length(frames1) ~= length(frames2)
    break
end

frames = frames1;
t = t(frames);
dt = t(2) - t(1);

x_head1(:,:) = data.kine.FLYBODY1.data.coords(1,1,frames);
y_head1(:,:) = data.kine.FLYBODY1.data.coords(1,2,frames);
z_head1(:,:) = data.kine.FLYBODY1.data.coords(1,3,frames);

x_tail1(:,:) = data.kine.FLYBODY1.data.coords(2,1,frames);
y_tail1(:,:) = data.kine.FLYBODY1.data.coords(2,2,frames);
z_tail1(:,:) = data.kine.FLYBODY1.data.coords(2,3,frames);

x_head2(:,:) = data.kine.FLYBODY2.data.coords(1,1,frames);
y_head2(:,:) = data.kine.FLYBODY2.data.coords(1,2,frames);
z_head2(:,:) = data.kine.FLYBODY2.data.coords(1,3,frames);

x_tail2(:,:) = data.kine.FLYBODY2.data.coords(2,1,frames);
y_tail2(:,:) = data.kine.FLYBODY2.data.coords(2,2,frames);
z_tail2(:,:) = data.kine.FLYBODY2.data.coords(2,3,frames);

d = sqrt( (x_head1 - x_head2).^2 + (y_head1 - y_head2).^2 + (z_head1 - z_head2).^2 );
teta = 2* atand(b/2./d);










