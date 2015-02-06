z(:,:) = pathDB.pos(:,:,3);
z_low = [];
z_mid = [];
z_high = [];
for i=1:size(z,2)
    z_temp = z(:,i);
    z_temp = z_temp(isnan(z_temp)==0);
    z_temp = z_temp(end);
    
    if settings.expansion.VerPos(i) == -1
        z_low(length(z_low)+1) = z_temp;
    elseif settings.expansion.VerPos(i) == 0
        z_mid(length(z_mid)+1) = z_temp;
    elseif settings.expansion.VerPos(i) == 1
        z_high(length(z_high)+1) = z_temp;
    end
    
end


errorbar(-53,mean(z_low),std(z_low),'xb')
hold on
errorbar(0,mean(z_mid),std(z_mid),'xr')
errorbar(53,mean(z_high),std(z_high),'xg')
legend('low','mid','high')

plot(-53,z_low,'ok')
plot(0,z_mid,'ok')
plot(53,z_high,'ok')

errorbar(-53,mean(z_low),std(z_low),'xb')
errorbar(0,mean(z_mid),std(z_mid),'xr')
errorbar(53,mean(z_high),std(z_high),'xg')

grid on
axis([-60 60 0 0.05])
xlabel('stimulus angle')
ylabel('escape height')

