
for i = 1:size(pathDB.pos,2)
    if toplot == 1

        t = pathDB.t;
        
        steady_frames = pathDB.steady_frames(:,i);

        stim_angle_vel = pathDB.stim_angle_vel(:,i);
        stim_angle_accel = pathDB.stim_angle_accel(:,i);

        slip = pathDB.slip(:,i);
        pitch = pathDB.pitch(:,i);
        roll = pathDB.roll(:,i);

        V = pathDB.V(:,i);
        An_hor = pathDB.An_hor(:,i);
        At_hor = pathDB.At_hor(:,i);
        
        A = pathDB.A(:,i);
        A_hor = pathDB.A_hor(:,i);
        A_ver = pathDB.A_ver(:,i);

        plot_flighttracks_clusters_separate_n_CLIP
        cd(figdir)
        saveas(gcf,['flightpath',num2str(i),'_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.fig'])
        saveas(gcf,['flightpath',num2str(i),'_',num2str(settings.seq(i,1)),'_seq',num2str(settings.seq(i,2)),'.png'])
        cd ..
%         pause
    end
end

    
