% make borf plots F&M timelines

steps = [0:1/Nsteps:1];
steps = steps(1:end-1);

plot([-80 -80+Fref*Fscale],[-10 -10],'-k','linewidth',1)
hold on
plot([-10 10],[0 0],'-','linewidth',.5,'color',[.75 .75 .75])
plot([0 0],[-10 10],'-','linewidth',.5,'color',[.75 .75 .75])
axis equal
axis([-90 90 -45 45]) 
axis off

%% damaged wing
F_color = [1 0 0];
cord_color = [0 0 0];
    
    Fx_now = Fx_damagedwing_all(:,end);
    Fz_now = Fz_damagedwing_all(:,end);
    
    stroke_now = stroke_damaged_all(:,end);
    rot_now = -rot_damaged_all(:,end)+90;
    dev_now = dev_damaged_all(:,end);
    
    t_now = [0:1/(length(Fx_now)-1):1];
    
    plot(stroke_now,dev_now,'-','linewidth',.5,'color',cord_color)
    hold on
    
    for j = 1:length(steps)
        
        N = round(steps(j)*(length(Fx_now)-1) +1);
        
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_damaged/2 * sind(rot_now(N));
        dy = c_damaged/2 * cosd(rot_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',3,'markerfacecolor',cord_color,'color',cord_color)
    end

%% intact wing
F_color = [0 0 1];
cord_color = [.5 .5 .5];
    
    Fx_now = Fx_intactwing_all(:,end);
    Fz_now = Fz_intactwing_all(:,end);
    
    stroke_now = stroke_intact_all(:,end);
    rot_now = -rot_intact_all(:,end)+90;
    dev_now = dev_intact_all(:,end);
    
    t_now = [0:1/(length(Fx_now)-1):1];

    plot(stroke_now,dev_now,'-','linewidth',.5,'color',cord_color)
    hold on
    
    for j = 1:length(steps)
        
        N = round(steps(j)*(length(Fx_now)-1) +1);
        
        x_wing = stroke_now(N);
        y_wing = dev_now(N);
        dx = c_intact/2 * sind(rot_now(N));
        dy = c_intact/2 * cosd(rot_now(N));
        
        x_LE = x_wing + dx;
        x_TE = x_wing - dx;
        y_LE = y_wing + dy;
        y_TE = y_wing - dy;
        
        Fx_tip = x_wing - Fscale*Fx_now(N);
        Fy_tip = y_wing - Fscale*Fz_now(N);
        
        plot([x_wing Fx_tip],[y_wing Fy_tip],'-','linewidth',1,'color',F_color)
        plot([x_LE x_TE],[y_LE y_TE],'-','linewidth',2,'color',cord_color)
        plot([x_LE],[y_LE],'o','markersize',3,'markerfacecolor',cord_color,'color',cord_color)
    end
