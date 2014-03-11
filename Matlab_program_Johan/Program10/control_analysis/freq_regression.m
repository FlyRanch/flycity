function freq_regression( non_dim_mat )


    freq_Fz_fit = polyfit(non_dim_mat(:,1),sqrt(non_dim_mat(:,3).^2+non_dim_mat(:,4).^2+non_dim_mat(:,5).^2)+sqrt(non_dim_mat(:,6).^2+non_dim_mat(:,7).^2+non_dim_mat(:,8).^2),1);
    
    downup_Fz_fit = polyfit(non_dim_mat(:,2),sqrt(non_dim_mat(:,3).^2+non_dim_mat(:,4).^2+non_dim_mat(:,5).^2)+sqrt(non_dim_mat(:,6).^2+non_dim_mat(:,7).^2+non_dim_mat(:,8).^2),1);
    
    freq_downup = polyfit(non_dim_mat(:,1),non_dim_mat(:,2),2);

    figure()
    hold on
    plot(non_dim_mat(:,1),sqrt(non_dim_mat(:,3).^2+non_dim_mat(:,4).^2+non_dim_mat(:,5).^2)+sqrt(non_dim_mat(:,6).^2+non_dim_mat(:,7).^2+non_dim_mat(:,8).^2),'o')
    plot(sort(non_dim_mat(:,1)),polyval(freq_Fz_fit,sort(non_dim_mat(:,1))),'r')
    hold off
    
    figure()
    hold on
    plot(non_dim_mat(:,2),sqrt(non_dim_mat(:,3).^2+non_dim_mat(:,4).^2+non_dim_mat(:,5).^2)+sqrt(non_dim_mat(:,6).^2+non_dim_mat(:,7).^2+non_dim_mat(:,8).^2),'o')
    plot(sort(non_dim_mat(:,2)),polyval(downup_Fz_fit,sort(non_dim_mat(:,2))),'r')
    hold off
    
    figure()
    hold on
    plot(non_dim_mat(:,1),non_dim_mat(:,2),'o')
    plot(sort(non_dim_mat(:,1)),polyval(freq_downup,sort(non_dim_mat(:,1))),'r')
    hold off

end

