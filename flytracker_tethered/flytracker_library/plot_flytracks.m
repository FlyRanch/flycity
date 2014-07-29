nframes = 400;
data = zeros(15,200);
for i = 1:nframes
    flstrng = sprintf('/home/matt/Dropbox/WORK/flytracker_thad/test_flytracks/fly%03d', i+6);
    tmp = load(flstrng);
    data(:,i) = tmp.xh;
end

rwing = data(12:12+3,:);
lwing = data(8:8+3,:);

lchord = [-2,-20,0;2,-20,0]';
rchord = [-2,20,0;2,20,0]';

ltip = [0,-20,0]';
rtip = [0,20,0]';

%[lroll,lpitch,lyaw] = qbody2angles_THAD(lwing);
%[rroll,rpitch,ryaw] = qbody2angles_THAD(rwing);

plot3(lchord(1,:),lchord(2,:),lchord(3,:),'g');
hold on;
plot3(rchord(1,:),rchord(2,:),rchord(3,:),'g');

for i = 1:nframes;
    q1 = rwing(:,i);
    rmat = quat2matNEW(q1);
    p_chord = rmat*rchord; %+ [0,5,0;0,5,0]';
    hold on;
    tip = rmat*rtip;
    if [0,1,0]*tip > 5
        plot3(p_chord(1,:),p_chord(2,:),p_chord(3,:),'b');
    else
        plot3(p_chord(1,:),p_chord(2,:),p_chord(3,:),'r');
    end
end

for i = 1:nframes;
    q1 = lwing(:,i);
    rmat = quat2matNEW(q1);
    p_chord = rmat*lchord; % + [0,-5,0;0,-5,0]';
    tip = rmat*ltip;
    if [0,1,0]*tip < -5
        plot3(p_chord(1,:),p_chord(2,:),p_chord(3,:),'b');
    else
        plot3(p_chord(1,:),p_chord(2,:),p_chord(3,:),'r');
    end
end

axis equal

rmats = nan(3,3,nframes);
tips = nan(3,nframes);

for i = 1:nframes
    rmats(:,:,i) = quat2matNEW(rwing(:,i));
    tips(:,i) = rmats(:,:,i)*rtip;
end

figure()
ea = SpinConv('DCMtoEA213',rmats,'1',1);
a = rad2deg(unwrap(deg2rad(ea(:,1))));
b = rad2deg(unwrap(deg2rad(ea(:,2))));
c = rad2deg(unwrap(deg2rad(ea(:,3))));

idx = tips(2,:)<0;

subplot(3,1,1);
plot(a,'-o');
a(idx) = -90;
hold on;
plot(a);
hold on; plot((tips(2,:)<0)*100);
subplot(3,1,2);
plot(b,'-o');
hold on; plot((tips(2,:)<0)*100);
subplot(3,1,3);
plot(c,'-o');
hold on; plot((tips(2,:)<0)*100);
plot(ones(100)*180);

c(idx) = -90;
hold on; plot(c,'r');

ea_prime = [a b c];
rmats_prime = SpinConv('EA213toDCM',ea_prime,'1',1);
for i = 1:nframes;
    tips_prime(:,i) = rmats_prime(:,:,i)*rtip;
end

%plot3(tips_prime(1,:),tips_prime(2,:),tips_prime(3,:),'r')
%%%%%%%%%%%%%
% figure()
% 
% for i = 1:100
%     rmats(:,:,i) = quat2matNEW(lwing(:,i));
%     tips(:,i) = rmats(:,:,i)*ltip;
% end
% 
% ea = SpinCalc('DCMtoEA213',rmats,'1',1);
% a = rad2deg(unwrap(deg2rad(ea(:,1))));
% b = rad2deg(unwrap(deg2rad(ea(:,2))));
% c = rad2deg(unwrap(deg2rad(ea(:,3))));
% 
% idx = tips(2,:)>0;
% 
% subplot(3,1,1);
% plot(a,'-o');
% a(idx) = -90;
% hold on; plot(a,'r')
% hold on; plot((tips(2,:)<0)*100);
% subplot(3,1,2);
% plot(b,'-o');
% hold on; plot((tips(2,:)<0)*100);
% subplot(3,1,3);
% plot(c,'-o');
% hold on; plot((tips(2,:)<0)*100);
% c(idx) = 90;
% hold on; plot(c,'r');
% 
% ea_prime = [a b c];
% rmats_prime = SpinCalc('EA213toDCM',ea_prime,'1',1);
% for i = 1:100;
%     tips_prime(:,i) = rmats_prime(:,:,i)*ltip;
% end
