
clear
clc
close all

load('cmap_grey2darkred.mat')
load('expansion_steps_Rlimited_exactpixval.mat')

clear t_new teta_new
t = [exp_steps(1,4)-exp_steps(2,4);exp_steps(:,4)];
teta = [0;exp_steps(:,3)];

for i=2:length(t);
    j=2*i-1;
    t_new(j-1:j,1)=t(i);
    teta_new(j:j+1,1)=teta(i);
end
t_new(j+1,1)=t(end);

figure(1)
plot(t_new,teta_new,'.-','color',cmap(1,:))

figure(2)
Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'-o');
for i=1:length(Y)
    set(Y(i),'color',cmap(1,:))
end 


load('expansion_steps_DOUBLESPEED_looped_exactpixval.mat')

clear t_new teta_new
t = [exp_steps(1,4)-exp_steps(2,4);exp_steps(:,4)];
teta = [0;exp_steps(:,3)];

for i=2:length(t);
    j=2*i-1;
    t_new(j-1:j,1)=t(i);
    teta_new(j:j+1,1)=teta(i);
end
t_new(j+1,1)=t(end);

figure(1)
hold on
plot(t_new,teta_new,'.-','color',cmap(end,:))

figure(2)
hold on
Y = polar(deg2rad(teta_new(1:14)/2),20*t_new(1:14)+10,'-o');
for i=1:length(Y)
    set(Y(i),'color',cmap(end,:))
end 

figure(1)
xlabel('time')
ylabel('optical angle')
set(gca,'xlim',[-0.02 .2],'ylim',[0 180])
set(gca,'XTick',[0:.05:.2])
set(gca,'YTick',[0:45:180])
saveas(gca,'expansion_angle_vs_time.fig')
saveas(gca,'expansion_angle_vs_time.jpg')

