% wrap data
if isempty(t_steady)==0
    t_steady_wrap = t_steady;
    stroke_wb_L_steady_wrap = stroke_wb_L_steady;
    stroke_wb_R_steady_wrap = stroke_wb_R_steady;
    pitch_wb_L_steady_wrap = pitch_wb_L_steady;
    pitch_wb_R_steady_wrap = pitch_wb_R_steady;
    dev_wb_L_steady_wrap = dev_wb_L_steady;
    dev_wb_R_steady_wrap = dev_wb_R_steady;
    aoa_wb_L_steady_wrap = aoa_wb_L_steady;
    aoa_wb_R_steady_wrap = aoa_wb_R_steady;
    U_wb_L_steady_wrap = U_wb_L_steady;
    U_wb_R_steady_wrap = U_wb_R_steady;

    Dstroke_wb_steady_wrap = Dstroke_wb_steady;
    Dpitch_wb_steady_wrap = Dpitch_wb_steady;
    Ddev_wb_steady_wrap = Ddev_wb_steady;
    Daoa_wb_steady_wrap = Daoa_wb_steady;
    DU_wb_steady_wrap = DU_wb_steady;

    t_steady_wrap(end+1:end+length(t_steady)) = t_steady -1;
    t_steady_wrap(end+1:end+length(t_steady)) = t_steady +1;
    
    stroke_wb_L_steady_wrap(end+1:end+length(stroke_wb_L_steady)) = stroke_wb_L_steady -1;
    stroke_wb_L_steady_wrap(end+1:end+length(stroke_wb_L_steady)) = stroke_wb_L_steady +1;
    pitch_wb_L_steady_wrap(end+1:end+length(pitch_wb_L_steady)) = pitch_wb_L_steady -1;
    pitch_wb_L_steady_wrap(end+1:end+length(pitch_wb_L_steady)) = pitch_wb_L_steady +1;
    dev_wb_L_steady_wrap(end+1:end+length(dev_wb_L_steady)) = dev_wb_L_steady -1;
    dev_wb_L_steady_wrap(end+1:end+length(dev_wb_L_steady)) = dev_wb_L_steady +1;
    aoa_wb_L_steady_wrap(end+1:end+length(aoa_wb_L_steady)) = aoa_wb_L_steady -1;
    aoa_wb_L_steady_wrap(end+1:end+length(aoa_wb_L_steady)) = aoa_wb_L_steady +1;
    U_wb_L_steady_wrap(end+1:end+length(U_wb_L_steady)) = U_wb_L_steady -1;
    U_wb_L_steady_wrap(end+1:end+length(U_wb_L_steady)) = U_wb_L_steady +1;
    
    stroke_wb_R_steady_wrap(end+1:end+length(stroke_wb_R_steady)) = stroke_wb_R_steady -1;
    stroke_wb_R_steady_wrap(end+1:end+length(stroke_wb_R_steady)) = stroke_wb_R_steady +1;
    pitch_wb_R_steady_wrap(end+1:end+length(pitch_wb_R_steady)) = pitch_wb_R_steady -1;
    pitch_wb_R_steady_wrap(end+1:end+length(pitch_wb_R_steady)) = pitch_wb_R_steady +1;
    dev_wb_R_steady_wrap(end+1:end+length(dev_wb_R_steady)) = dev_wb_R_steady -1;
    dev_wb_R_steady_wrap(end+1:end+length(dev_wb_R_steady)) = dev_wb_R_steady +1;
    aoa_wb_R_steady_wrap(end+1:end+length(aoa_wb_R_steady)) = aoa_wb_R_steady -1;
    aoa_wb_R_steady_wrap(end+1:end+length(aoa_wb_R_steady)) = aoa_wb_R_steady +1;
    U_wb_R_steady_wrap(end+1:end+length(U_wb_R_steady)) = U_wb_R_steady -1;
    U_wb_R_steady_wrap(end+1:end+length(U_wb_R_steady)) = U_wb_R_steady +1;

    Dstroke_wb_steady_wrap(end+1:end+length(Dstroke_wb_steady)) = Dstroke_wb_steady -1;
    Dstroke_wb_steady_wrap(end+1:end+length(Dstroke_wb_steady)) = Dstroke_wb_steady +1;
    Dpitch_wb_steady_wrap(end+1:end+length(Dpitch_wb_steady)) = Dpitch_wb_steady -1;
    Dpitch_wb_steady_wrap(end+1:end+length(Dpitch_wb_steady)) = Dpitch_wb_steady +1;
    Ddev_wb_steady_wrap(end+1:end+length(Ddev_wb_steady)) = Ddev_wb_steady -1;
    Ddev_wb_steady_wrap(end+1:end+length(Ddev_wb_steady)) = Ddev_wb_steady +1;
    Daoa_wb_steady_wrap(end+1:end+length(Daoa_wb_steady)) = Daoa_wb_steady -1;
    Daoa_wb_steady_wrap(end+1:end+length(Daoa_wb_steady)) = Daoa_wb_steady +1;
    DU_wb_steady_wrap(end+1:end+length(DU_wb_steady)) = DU_wb_steady -1;
    DU_wb_steady_wrap(end+1:end+length(DU_wb_steady)) = DU_wb_steady +1;
    
    stroke_wb_L_steady_func = csaps(t_steady_wrap,stroke_wb_L_steady_wrap);
    pitch_wb_L_steady_func = csaps(t_steady_wrap,pitch_wb_L_steady_wrap);
    dev_wb_L_steady_func = csaps(t_steady_wrap,dev_wb_L_steady_wrap);
    aoa_wb_L_steady_func = csaps(t_steady_wrap,aoa_wb_L_steady_wrap);
    U_wb_L_steady_func = csaps(t_steady_wrap,U_wb_L_steady_wrap);
    
    stroke_wb_R_steady_func = csaps(t_steady_wrap,stroke_wb_R_steady_wrap);
    pitch_wb_R_steady_func = csaps(t_steady_wrap,pitch_wb_R_steady_wrap);
    dev_wb_R_steady_func = csaps(t_steady_wrap,dev_wb_R_steady_wrap);
    aoa_wb_R_steady_func = csaps(t_steady_wrap,aoa_wb_R_steady_wrap);
    U_wb_R_steady_func = csaps(t_steady_wrap,U_wb_R_steady_wrap);
    
    Dstroke_wb_steady_func = csaps(t_steady_wrap,Dstroke_wb_steady_wrap);
    Dpitch_wb_steady_func = csaps(t_steady_wrap,Dpitch_wb_steady_wrap);
    Ddev_wb_steady_func = csaps(t_steady_wrap,Ddev_wb_steady_wrap);
    Daoa_wb_steady_func = csaps(t_steady_wrap,Daoa_wb_steady_wrap);
    DU_wb_steady_func = csaps(t_steady_wrap,DU_wb_steady_wrap);
end
    
if isempty(t_hi)==0
    t_hi_wrap = t_hi;
    stroke_wb_L_hi_wrap = stroke_wb_L_hi;
    stroke_wb_R_hi_wrap = stroke_wb_R_hi;
    pitch_wb_L_hi_wrap = pitch_wb_L_hi;
    pitch_wb_R_hi_wrap = pitch_wb_R_hi;
    dev_wb_L_hi_wrap = dev_wb_L_hi;
    dev_wb_R_hi_wrap = dev_wb_R_hi;
    aoa_wb_L_hi_wrap = aoa_wb_L_hi;
    aoa_wb_R_hi_wrap = aoa_wb_R_hi;
    U_wb_L_hi_wrap = U_wb_L_hi;
    U_wb_R_hi_wrap = U_wb_R_hi;

    Dstroke_wb_hi_wrap = Dstroke_wb_hi;
    Dpitch_wb_hi_wrap = Dpitch_wb_hi;
    Ddev_wb_hi_wrap = Ddev_wb_hi;
    Daoa_wb_hi_wrap = Daoa_wb_hi;
    DU_wb_hi_wrap = DU_wb_hi;
    
    t_hi_wrap(end+1:end+length(t_hi)) = t_hi -1;
    t_hi_wrap(end+1:end+length(t_hi)) = t_hi +1;
    stroke_wb_L_hi_wrap(end+1:end+length(stroke_wb_L_hi)) = stroke_wb_L_hi -1;
    stroke_wb_L_hi_wrap(end+1:end+length(stroke_wb_L_hi)) = stroke_wb_L_hi +1;
    pitch_wb_L_hi_wrap(end+1:end+length(pitch_wb_L_hi)) = pitch_wb_L_hi -1;
    pitch_wb_L_hi_wrap(end+1:end+length(pitch_wb_L_hi)) = pitch_wb_L_hi +1;
    dev_wb_L_hi_wrap(end+1:end+length(dev_wb_L_hi)) = dev_wb_L_hi -1;
    dev_wb_L_hi_wrap(end+1:end+length(dev_wb_L_hi)) = dev_wb_L_hi +1;
    aoa_wb_L_hi_wrap(end+1:end+length(aoa_wb_L_hi)) = aoa_wb_L_hi -1;
    aoa_wb_L_hi_wrap(end+1:end+length(aoa_wb_L_hi)) = aoa_wb_L_hi +1;
    U_wb_L_hi_wrap(end+1:end+length(U_wb_L_hi)) = U_wb_L_hi -1;
    U_wb_L_hi_wrap(end+1:end+length(U_wb_L_hi)) = U_wb_L_hi +1;
    
    stroke_wb_R_hi_wrap(end+1:end+length(stroke_wb_R_hi)) = stroke_wb_R_hi -1;
    stroke_wb_R_hi_wrap(end+1:end+length(stroke_wb_R_hi)) = stroke_wb_R_hi +1;
    pitch_wb_R_hi_wrap(end+1:end+length(pitch_wb_R_hi)) = pitch_wb_R_hi -1;
    pitch_wb_R_hi_wrap(end+1:end+length(pitch_wb_R_hi)) = pitch_wb_R_hi +1;
    dev_wb_R_hi_wrap(end+1:end+length(dev_wb_R_hi)) = dev_wb_R_hi -1;
    dev_wb_R_hi_wrap(end+1:end+length(dev_wb_R_hi)) = dev_wb_R_hi +1;
    aoa_wb_R_hi_wrap(end+1:end+length(aoa_wb_R_hi)) = aoa_wb_R_hi -1;
    aoa_wb_R_hi_wrap(end+1:end+length(aoa_wb_R_hi)) = aoa_wb_R_hi +1;
    U_wb_R_hi_wrap(end+1:end+length(U_wb_R_hi)) = U_wb_R_hi -1;
    U_wb_R_hi_wrap(end+1:end+length(U_wb_R_hi)) = U_wb_R_hi +1;

    Dstroke_wb_hi_wrap(end+1:end+length(Dstroke_wb_hi)) = Dstroke_wb_hi -1;
    Dstroke_wb_hi_wrap(end+1:end+length(Dstroke_wb_hi)) = Dstroke_wb_hi +1;
    Dpitch_wb_hi_wrap(end+1:end+length(Dpitch_wb_hi)) = Dpitch_wb_hi -1;
    Dpitch_wb_hi_wrap(end+1:end+length(Dpitch_wb_hi)) = Dpitch_wb_hi +1;
    Ddev_wb_hi_wrap(end+1:end+length(Ddev_wb_hi)) = Ddev_wb_hi -1;
    Ddev_wb_hi_wrap(end+1:end+length(Ddev_wb_hi)) = Ddev_wb_hi +1;
    Daoa_wb_hi_wrap(end+1:end+length(Daoa_wb_hi)) = Daoa_wb_hi -1;
    Daoa_wb_hi_wrap(end+1:end+length(Daoa_wb_hi)) = Daoa_wb_hi +1;
    DU_wb_hi_wrap(end+1:end+length(DU_wb_hi)) = DU_wb_hi -1;
    DU_wb_hi_wrap(end+1:end+length(DU_wb_hi)) = DU_wb_hi +1;
    
    stroke_wb_L_hi_func = csaps(t_hi_wrap,stroke_wb_L_hi_wrap);
    pitch_wb_L_hi_func = csaps(t_hi_wrap,pitch_wb_L_hi_wrap);
    dev_wb_L_hi_func = csaps(t_hi_wrap,dev_wb_L_hi_wrap);
    aoa_wb_L_hi_func = csaps(t_hi_wrap,aoa_wb_L_hi_wrap);
    U_wb_L_hi_func = csaps(t_hi_wrap,U_wb_L_hi_wrap);
    
    stroke_wb_R_hi_func = csaps(t_hi_wrap,stroke_wb_R_hi_wrap);
    pitch_wb_R_hi_func = csaps(t_hi_wrap,pitch_wb_R_hi_wrap);
    dev_wb_R_hi_func = csaps(t_hi_wrap,dev_wb_R_hi_wrap);
    aoa_wb_R_hi_func = csaps(t_hi_wrap,aoa_wb_R_hi_wrap);
    U_wb_R_hi_func = csaps(t_hi_wrap,U_wb_R_hi_wrap);
        
    Dstroke_wb_hi_func = csaps(t_hi_wrap,Dstroke_wb_hi_wrap);
    Dpitch_wb_hi_func = csaps(t_hi_wrap,Dpitch_wb_hi_wrap);
    Ddev_wb_hi_func = csaps(t_hi_wrap,Ddev_wb_hi_wrap);
    Daoa_wb_hi_func = csaps(t_hi_wrap,Daoa_wb_hi_wrap);
    DU_wb_hi_func = csaps(t_hi_wrap,DU_wb_hi_wrap);

end

if isempty(t_low)==0
    t_low_wrap = t_low;
    stroke_wb_L_low_wrap = stroke_wb_L_low;
    stroke_wb_R_low_wrap = stroke_wb_R_low;
    pitch_wb_L_low_wrap = pitch_wb_L_low;
    pitch_wb_R_low_wrap = pitch_wb_R_low;
    dev_wb_L_low_wrap = dev_wb_L_low;
    dev_wb_R_low_wrap = dev_wb_R_low;
    aoa_wb_L_low_wrap = aoa_wb_L_low;
    aoa_wb_R_low_wrap = aoa_wb_R_low;
    U_wb_L_low_wrap = U_wb_L_low;
    U_wb_R_low_wrap = U_wb_R_low;

    Dstroke_wb_low_wrap = Dstroke_wb_low;
    Dpitch_wb_low_wrap = Dpitch_wb_low;
    Ddev_wb_low_wrap = Ddev_wb_low;
    Daoa_wb_low_wrap = Daoa_wb_low;
    DU_wb_low_wrap = DU_wb_low;

    t_low_wrap(end+1:end+length(t_low)) = t_low -1;
    t_low_wrap(end+1:end+length(t_low)) = t_low +1;
    stroke_wb_L_low_wrap(end+1:end+length(stroke_wb_L_low)) = stroke_wb_L_low -1;
    stroke_wb_L_low_wrap(end+1:end+length(stroke_wb_L_low)) = stroke_wb_L_low +1;
    pitch_wb_L_low_wrap(end+1:end+length(pitch_wb_L_low)) = pitch_wb_L_low -1;
    pitch_wb_L_low_wrap(end+1:end+length(pitch_wb_L_low)) = pitch_wb_L_low +1;
    dev_wb_L_low_wrap(end+1:end+length(dev_wb_L_low)) = dev_wb_L_low -1;
    dev_wb_L_low_wrap(end+1:end+length(dev_wb_L_low)) = dev_wb_L_low +1;
    aoa_wb_L_low_wrap(end+1:end+length(aoa_wb_L_low)) = aoa_wb_L_low -1;
    aoa_wb_L_low_wrap(end+1:end+length(aoa_wb_L_low)) = aoa_wb_L_low +1;
    U_wb_L_low_wrap(end+1:end+length(U_wb_L_low)) = U_wb_L_low -1;
    U_wb_L_low_wrap(end+1:end+length(U_wb_L_low)) = U_wb_L_low +1;
    
    stroke_wb_R_low_wrap(end+1:end+length(stroke_wb_R_low)) = stroke_wb_R_low -1;
    stroke_wb_R_low_wrap(end+1:end+length(stroke_wb_R_low)) = stroke_wb_R_low +1;
    pitch_wb_R_low_wrap(end+1:end+length(pitch_wb_R_low)) = pitch_wb_R_low -1;
    pitch_wb_R_low_wrap(end+1:end+length(pitch_wb_R_low)) = pitch_wb_R_low +1;
    dev_wb_R_low_wrap(end+1:end+length(dev_wb_R_low)) = dev_wb_R_low -1;
    dev_wb_R_low_wrap(end+1:end+length(dev_wb_R_low)) = dev_wb_R_low +1;
    aoa_wb_R_low_wrap(end+1:end+length(aoa_wb_R_low)) = aoa_wb_R_low -1;
    aoa_wb_R_low_wrap(end+1:end+length(aoa_wb_R_low)) = aoa_wb_R_low +1;
    U_wb_R_low_wrap(end+1:end+length(U_wb_R_low)) = U_wb_R_low -1;
    U_wb_R_low_wrap(end+1:end+length(U_wb_R_low)) = U_wb_R_low +1;

    Dstroke_wb_low_wrap(end+1:end+length(Dstroke_wb_low)) = Dstroke_wb_low -1;
    Dstroke_wb_low_wrap(end+1:end+length(Dstroke_wb_low)) = Dstroke_wb_low +1;
    Dpitch_wb_low_wrap(end+1:end+length(Dpitch_wb_low)) = Dpitch_wb_low -1;
    Dpitch_wb_low_wrap(end+1:end+length(Dpitch_wb_low)) = Dpitch_wb_low +1;
    Ddev_wb_low_wrap(end+1:end+length(Ddev_wb_low)) = Ddev_wb_low -1;
    Ddev_wb_low_wrap(end+1:end+length(Ddev_wb_low)) = Ddev_wb_low +1;
    Daoa_wb_low_wrap(end+1:end+length(Daoa_wb_low)) = Daoa_wb_low -1;
    Daoa_wb_low_wrap(end+1:end+length(Daoa_wb_low)) = Daoa_wb_low +1;
    DU_wb_low_wrap(end+1:end+length(DU_wb_low)) = DU_wb_low -1;
    DU_wb_low_wrap(end+1:end+length(DU_wb_low)) = DU_wb_low +1;
    
    stroke_wb_L_low_func = csaps(t_low_wrap,stroke_wb_L_low_wrap);
    pitch_wb_L_low_func = csaps(t_low_wrap,pitch_wb_L_low_wrap);
    dev_wb_L_low_func = csaps(t_low_wrap,dev_wb_L_low_wrap);
    aoa_wb_L_low_func = csaps(t_low_wrap,aoa_wb_L_low_wrap);
    U_wb_L_low_func = csaps(t_low_wrap,U_wb_L_low_wrap);
    
    stroke_wb_R_low_func = csaps(t_low_wrap,stroke_wb_R_low_wrap);
    pitch_wb_R_low_func = csaps(t_low_wrap,pitch_wb_R_low_wrap);
    dev_wb_R_low_func = csaps(t_low_wrap,dev_wb_R_low_wrap);
    aoa_wb_R_low_func = csaps(t_low_wrap,aoa_wb_R_low_wrap);
    U_wb_R_low_func = csaps(t_low_wrap,U_wb_R_low_wrap);
        
    Dstroke_wb_low_func = csaps(t_low_wrap,Dstroke_wb_low_wrap);
    Dpitch_wb_low_func = csaps(t_low_wrap,Dpitch_wb_low_wrap);
    Ddev_wb_low_func = csaps(t_low_wrap,Ddev_wb_low_wrap);
    Daoa_wb_low_func = csaps(t_low_wrap,Daoa_wb_low_wrap);
    DU_wb_low_func = csaps(t_low_wrap,DU_wb_low_wrap);

end
