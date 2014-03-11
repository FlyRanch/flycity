function Flight_trajectory_plot( settings, pathDB )

    
    % Create the flight trajectory plot
    
%     figure(1)
%     hFig = figure(1);
%     set(gcf,'PaperPositionMode','auto');
%     set(hFig,'Position',[0 0 1000 1000]);
%     flight_traject( settings, pathDB, 71 )
%     title('Roll maneuver, seq 71')
%     saveas(hFig,[char(settings.plot_loc) '/Escape_maneuvers/Flight_traject_71'],'fig')
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 1000]);
    flight_traject( settings, pathDB, 58 )
    title('Roll maneuver, seq 58')
%     saveas(hFig,[char(settings.plot_loc) '/Escape_maneuvers/Flight_traject_71'],'fig')

    pause
    
    close all
        
%     % Plot wingkinematics related to the wingbeats:
%     
%     for i = 1:pathDB.wingbeats.nr_of_wb(71)
%         figure(i)
%         hFig = figure(i);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 600 600]);
%         wingkin_plot_3D( settings, pathDB, 71, i)
%         [~,h1] = suplabel(char(['Wing kinematics of wingbeat ' num2str(i)]), 't');
%         set(h1,'FontSize',12)
%         saveas(hFig,[char(settings.plot_loc) char(['/Escape_maneuvers/Wingkin_71_wb_' num2str(i)])],'fig')
%     end
%     
%     close all
%     
%     % Plot the forces related to the wingbeats:
%     
%     for i = 1:pathDB.wingbeats.nr_of_wb(71)
%         figure(i)
%         hFig = figure(i);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 600 600]);
%         Aero_force_plot_3D( settings, pathDB, 71, i )
%         [~,h1] = suplabel(char(['Aerodynamic forces at wingbeat ' num2str(i)]), 't');
%         set(h1,'FontSize',12)
%         saveas(hFig,[char(settings.plot_loc) char(['/Escape_maneuvers/Forceplot_71_wb_' num2str(i)])],'fig')
%     end

    for i = 1:pathDB.wingbeats.nr_of_wb(71)
        figure(i)
        hFig = figure(i);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 1200 800]);
        wing_kin_aero_plot( settings, pathDB, 71, i)
        [~,h1] = suplabel(char(['Wingbeat ' num2str(i)]), 't');
        set(h1,'FontSize',12)
%         saveas(hFig,[char(settings.plot_loc) char(['/Escape_maneuvers/wingkin_force_71_' num2str(i)])],'fig')
    end
    
    
    close all
    
       
    % Create the flight trajectory plot
    
    figure(1)
    hFig = figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(hFig,'Position',[0 0 1000 1000]);
    flight_traject( settings, pathDB, 99 )
    title('Pitch maneuver, seq 99')
    saveas(hFig,[char(settings.plot_loc) '/Escape_maneuvers/Flight_traject_99'],'fig')
    
    close all
    
%     % Plot wingkinematics related to the wingbeats:
%     
%     for i = 1:pathDB.wingbeats.nr_of_wb(99)
%         figure(i)
%         hFig = figure(i);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 600 600]);
%         wingkin_plot_3D( settings, pathDB, 99, i)
%         [~,h1] = suplabel(char(['Wing kinematics of wingbeat ' num2str(i)]), 't');
%         set(h1,'FontSize',12)
%         saveas(hFig,[char(settings.plot_loc) char(['/Escape_maneuvers/Wingkin_99_wb_' num2str(i)])],'fig')
%     end
%     
%     close all
%         
%     % Plot the forces related to the wingbeats:
%     
%     for i = 1:pathDB.wingbeats.nr_of_wb(99)
%         figure(i+1)
%         hFig = figure(i+1);
%         set(gcf,'PaperPositionMode','auto');
%         set(hFig,'Position',[0 0 600 600]);
%         Aero_force_plot_3D( settings, pathDB, 99, i )
%         [~,h1] = suplabel(char(['Aerodynamic forces at wingbeat ' num2str(i)]), 't');
%         set(h1,'FontSize',12)
%         saveas(hFig,[char(settings.plot_loc) char(['/Escape_maneuvers/Forceplot_99_wb_' num2str(i)])],'fig')
%     end

    for i = 1:pathDB.wingbeats.nr_of_wb(99)
        figure(i)
        hFig = figure(i);
        set(gcf,'PaperPositionMode','auto');
        set(hFig,'Position',[0 0 1200 800]);
        wing_kin_aero_plot( settings, pathDB, 99, i)
        [~,h1] = suplabel(char(['Wingbeat ' num2str(i)]), 't');
        set(h1,'FontSize',12)
        saveas(hFig,[char(settings.plot_loc) char(['/Escape_maneuvers/wingkin_force_99_' num2str(i)])],'fig')
    end
    
    close all

end

