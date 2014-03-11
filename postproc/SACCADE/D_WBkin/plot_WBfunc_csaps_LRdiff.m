if isempty(t_steady)==0
    plot_color = color_mid;
    
    subplot(3,5,1)
    yi = fnval(stroke_wb_L_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,6)
    yi = fnval(stroke_wb_R_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,11)
    yi = fnval(Dstroke_wb_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,2)
    yi = fnval(pitch_wb_L_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,7)
    yi = fnval(pitch_wb_R_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,12)
    yi = fnval(Dpitch_wb_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,3)
    yi = fnval(dev_wb_L_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,8)
    yi = fnval(dev_wb_R_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,13)
    yi = fnval(Ddev_wb_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,4)
    yi = fnval(aoa_wb_L_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,9)
    yi = fnval(aoa_wb_R_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,14)
    yi = fnval(Daoa_wb_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,5)
    yi = fnval(U_wb_L_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,10)
    yi = fnval(U_wb_R_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,15)
    yi = fnval(DU_wb_steady_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
end

if isempty(t_hi)==0
    plot_color = color_max;
    
    subplot(3,5,1)
    yi = fnval(stroke_wb_L_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,6)
    yi = fnval(stroke_wb_R_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,11)
    yi = fnval(Dstroke_wb_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,2)
    yi = fnval(pitch_wb_L_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,7)
    yi = fnval(pitch_wb_R_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,12)
    yi = fnval(Dpitch_wb_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,3)
    yi = fnval(dev_wb_L_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,8)
    yi = fnval(dev_wb_R_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,13)
    yi = fnval(Ddev_wb_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,4)
    yi = fnval(aoa_wb_L_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,9)
    yi = fnval(aoa_wb_R_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,14)
    yi = fnval(Daoa_wb_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,5)
    yi = fnval(U_wb_L_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,10)
    yi = fnval(U_wb_R_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,15)
    yi = fnval(DU_wb_hi_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
end

if isempty(t_low)==0
    plot_color = color_min;
    
    subplot(3,5,1)
    yi = fnval(stroke_wb_L_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,6)
    yi = fnval(stroke_wb_R_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,11)
    yi = fnval(Dstroke_wb_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,2)
    yi = fnval(pitch_wb_L_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,7)
    yi = fnval(pitch_wb_R_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,12)
    yi = fnval(Dpitch_wb_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,3)
    yi = fnval(dev_wb_L_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,8)
    yi = fnval(dev_wb_R_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,13)
    yi = fnval(Ddev_wb_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,4)
    yi = fnval(aoa_wb_L_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,9)
    yi = fnval(aoa_wb_R_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,14)
    yi = fnval(Daoa_wb_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    
    subplot(3,5,5)
    yi = fnval(U_wb_L_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,10)
    yi = fnval(U_wb_R_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
    subplot(3,5,15)
    yi = fnval(DU_wb_low_func,xi);
    plot(xi,yi,'-','color',plot_color,'linewidth',linewidth_meanWBs)
end

