% wrap data
if isempty(t_steady_bins(:))==0
    t_steady_wrap = t_steady_bins(:);
    stroke_wb_L_steady_wrap = stroke_wb_L_steady_bins(:);
    stroke_wb_R_steady_wrap = stroke_wb_R_steady_bins(:);
    pitch_wb_L_steady_wrap = pitch_wb_L_steady_bins(:);
    pitch_wb_R_steady_wrap = pitch_wb_R_steady_bins(:);
    dev_wb_L_steady_wrap = dev_wb_L_steady_bins(:);
    dev_wb_R_steady_wrap = dev_wb_R_steady_bins(:);
    aoa_wb_L_steady_wrap = aoa_wb_L_steady_bins(:);
    aoa_wb_R_steady_wrap = aoa_wb_R_steady_bins(:);
    U_wb_L_steady_wrap = U_wb_L_steady_bins(:);
    U_wb_R_steady_wrap = U_wb_R_steady_bins(:);

    Dstroke_wb_steady_wrap = Dstroke_wb_steady_bins(:);
    Dpitch_wb_steady_wrap = Dpitch_wb_steady_bins(:);
    Ddev_wb_steady_wrap = Ddev_wb_steady_bins(:);
    Daoa_wb_steady_wrap = Daoa_wb_steady_bins(:);
    DU_wb_steady_wrap = DU_wb_steady_bins(:);

    t_steady_wrap(end+1:end+length(t_steady_bins(:))) = t_steady_bins(:) -1;
    t_steady_wrap(end+1:end+length(t_steady_bins(:))) = t_steady_bins(:) +1;
    
    stroke_wb_L_steady_wrap(end+1:end+length(stroke_wb_L_steady_bins(:))) = stroke_wb_L_steady_bins(:) -1;
    stroke_wb_L_steady_wrap(end+1:end+length(stroke_wb_L_steady_bins(:))) = stroke_wb_L_steady_bins(:) +1;
    pitch_wb_L_steady_wrap(end+1:end+length(pitch_wb_L_steady_bins(:))) = pitch_wb_L_steady_bins(:) -1;
    pitch_wb_L_steady_wrap(end+1:end+length(pitch_wb_L_steady_bins(:))) = pitch_wb_L_steady_bins(:) +1;
    dev_wb_L_steady_wrap(end+1:end+length(dev_wb_L_steady_bins(:))) = dev_wb_L_steady_bins(:) -1;
    dev_wb_L_steady_wrap(end+1:end+length(dev_wb_L_steady_bins(:))) = dev_wb_L_steady_bins(:) +1;
    aoa_wb_L_steady_wrap(end+1:end+length(aoa_wb_L_steady_bins(:))) = aoa_wb_L_steady_bins(:) -1;
    aoa_wb_L_steady_wrap(end+1:end+length(aoa_wb_L_steady_bins(:))) = aoa_wb_L_steady_bins(:) +1;
    U_wb_L_steady_wrap(end+1:end+length(U_wb_L_steady_bins(:))) = U_wb_L_steady_bins(:) -1;
    U_wb_L_steady_wrap(end+1:end+length(U_wb_L_steady_bins(:))) = U_wb_L_steady_bins(:) +1;
    
    stroke_wb_R_steady_wrap(end+1:end+length(stroke_wb_R_steady_bins(:))) = stroke_wb_R_steady_bins(:) -1;
    stroke_wb_R_steady_wrap(end+1:end+length(stroke_wb_R_steady_bins(:))) = stroke_wb_R_steady_bins(:) +1;
    pitch_wb_R_steady_wrap(end+1:end+length(pitch_wb_R_steady_bins(:))) = pitch_wb_R_steady_bins(:) -1;
    pitch_wb_R_steady_wrap(end+1:end+length(pitch_wb_R_steady_bins(:))) = pitch_wb_R_steady_bins(:) +1;
    dev_wb_R_steady_wrap(end+1:end+length(dev_wb_R_steady_bins(:))) = dev_wb_R_steady_bins(:) -1;
    dev_wb_R_steady_wrap(end+1:end+length(dev_wb_R_steady_bins(:))) = dev_wb_R_steady_bins(:) +1;
    aoa_wb_R_steady_wrap(end+1:end+length(aoa_wb_R_steady_bins(:))) = aoa_wb_R_steady_bins(:) -1;
    aoa_wb_R_steady_wrap(end+1:end+length(aoa_wb_R_steady_bins(:))) = aoa_wb_R_steady_bins(:) +1;
    U_wb_R_steady_wrap(end+1:end+length(U_wb_R_steady_bins(:))) = U_wb_R_steady_bins(:) -1;
    U_wb_R_steady_wrap(end+1:end+length(U_wb_R_steady_bins(:))) = U_wb_R_steady_bins(:) +1;

    Dstroke_wb_steady_wrap(end+1:end+length(Dstroke_wb_steady_bins(:))) = Dstroke_wb_steady_bins(:) -1;
    Dstroke_wb_steady_wrap(end+1:end+length(Dstroke_wb_steady_bins(:))) = Dstroke_wb_steady_bins(:) +1;
    Dpitch_wb_steady_wrap(end+1:end+length(Dpitch_wb_steady_bins(:))) = Dpitch_wb_steady_bins(:) -1;
    Dpitch_wb_steady_wrap(end+1:end+length(Dpitch_wb_steady_bins(:))) = Dpitch_wb_steady_bins(:) +1;
    Ddev_wb_steady_wrap(end+1:end+length(Ddev_wb_steady_bins(:))) = Ddev_wb_steady_bins(:) -1;
    Ddev_wb_steady_wrap(end+1:end+length(Ddev_wb_steady_bins(:))) = Ddev_wb_steady_bins(:) +1;
    Daoa_wb_steady_wrap(end+1:end+length(Daoa_wb_steady_bins(:))) = Daoa_wb_steady_bins(:) -1;
    Daoa_wb_steady_wrap(end+1:end+length(Daoa_wb_steady_bins(:))) = Daoa_wb_steady_bins(:) +1;
    DU_wb_steady_wrap(end+1:end+length(DU_wb_steady_bins(:))) = DU_wb_steady_bins(:) -1;
    DU_wb_steady_wrap(end+1:end+length(DU_wb_steady_bins(:))) = DU_wb_steady_bins(:) +1;
    
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
    
    % L&R
    stroke_steady_func = csaps([t_steady_wrap;t_steady_wrap],[stroke_wb_L_steady_wrap;stroke_wb_R_steady_wrap]);
    pitch_steady_func = csaps([t_steady_wrap;t_steady_wrap],[pitch_wb_L_steady_wrap;pitch_wb_R_steady_wrap]);
    dev_steady_func = csaps([t_steady_wrap;t_steady_wrap],[dev_wb_L_steady_wrap;dev_wb_R_steady_wrap]);
    aoa_steady_func = csaps([t_steady_wrap;t_steady_wrap],[aoa_wb_L_steady_wrap;aoa_wb_R_steady_wrap]);
    U_steady_func = csaps([t_steady_wrap;t_steady_wrap],[U_wb_L_steady_wrap;U_wb_R_steady_wrap]);
end


if isempty(t_hi_bins(:))==0
    t_hi_wrap = t_hi_bins(:);
    stroke_wb_L_hi_wrap = stroke_wb_L_hi_bins(:);
    stroke_wb_R_hi_wrap = stroke_wb_R_hi_bins(:);
    pitch_wb_L_hi_wrap = pitch_wb_L_hi_bins(:);
    pitch_wb_R_hi_wrap = pitch_wb_R_hi_bins(:);
    dev_wb_L_hi_wrap = dev_wb_L_hi_bins(:);
    dev_wb_R_hi_wrap = dev_wb_R_hi_bins(:);
    aoa_wb_L_hi_wrap = aoa_wb_L_hi_bins(:);
    aoa_wb_R_hi_wrap = aoa_wb_R_hi_bins(:);
    U_wb_L_hi_wrap = U_wb_L_hi_bins(:);
    U_wb_R_hi_wrap = U_wb_R_hi_bins(:);

    Dstroke_wb_hi_wrap = Dstroke_wb_hi_bins(:);
    Dpitch_wb_hi_wrap = Dpitch_wb_hi_bins(:);
    Ddev_wb_hi_wrap = Ddev_wb_hi_bins(:);
    Daoa_wb_hi_wrap = Daoa_wb_hi_bins(:);
    DU_wb_hi_wrap = DU_wb_hi_bins(:);

    t_hi_wrap(end+1:end+length(t_hi_bins(:))) = t_hi_bins(:) -1;
    t_hi_wrap(end+1:end+length(t_hi_bins(:))) = t_hi_bins(:) +1;
    
    stroke_wb_L_hi_wrap(end+1:end+length(stroke_wb_L_hi_bins(:))) = stroke_wb_L_hi_bins(:) -1;
    stroke_wb_L_hi_wrap(end+1:end+length(stroke_wb_L_hi_bins(:))) = stroke_wb_L_hi_bins(:) +1;
    pitch_wb_L_hi_wrap(end+1:end+length(pitch_wb_L_hi_bins(:))) = pitch_wb_L_hi_bins(:) -1;
    pitch_wb_L_hi_wrap(end+1:end+length(pitch_wb_L_hi_bins(:))) = pitch_wb_L_hi_bins(:) +1;
    dev_wb_L_hi_wrap(end+1:end+length(dev_wb_L_hi_bins(:))) = dev_wb_L_hi_bins(:) -1;
    dev_wb_L_hi_wrap(end+1:end+length(dev_wb_L_hi_bins(:))) = dev_wb_L_hi_bins(:) +1;
    aoa_wb_L_hi_wrap(end+1:end+length(aoa_wb_L_hi_bins(:))) = aoa_wb_L_hi_bins(:) -1;
    aoa_wb_L_hi_wrap(end+1:end+length(aoa_wb_L_hi_bins(:))) = aoa_wb_L_hi_bins(:) +1;
    U_wb_L_hi_wrap(end+1:end+length(U_wb_L_hi_bins(:))) = U_wb_L_hi_bins(:) -1;
    U_wb_L_hi_wrap(end+1:end+length(U_wb_L_hi_bins(:))) = U_wb_L_hi_bins(:) +1;
    
    stroke_wb_R_hi_wrap(end+1:end+length(stroke_wb_R_hi_bins(:))) = stroke_wb_R_hi_bins(:) -1;
    stroke_wb_R_hi_wrap(end+1:end+length(stroke_wb_R_hi_bins(:))) = stroke_wb_R_hi_bins(:) +1;
    pitch_wb_R_hi_wrap(end+1:end+length(pitch_wb_R_hi_bins(:))) = pitch_wb_R_hi_bins(:) -1;
    pitch_wb_R_hi_wrap(end+1:end+length(pitch_wb_R_hi_bins(:))) = pitch_wb_R_hi_bins(:) +1;
    dev_wb_R_hi_wrap(end+1:end+length(dev_wb_R_hi_bins(:))) = dev_wb_R_hi_bins(:) -1;
    dev_wb_R_hi_wrap(end+1:end+length(dev_wb_R_hi_bins(:))) = dev_wb_R_hi_bins(:) +1;
    aoa_wb_R_hi_wrap(end+1:end+length(aoa_wb_R_hi_bins(:))) = aoa_wb_R_hi_bins(:) -1;
    aoa_wb_R_hi_wrap(end+1:end+length(aoa_wb_R_hi_bins(:))) = aoa_wb_R_hi_bins(:) +1;
    U_wb_R_hi_wrap(end+1:end+length(U_wb_R_hi_bins(:))) = U_wb_R_hi_bins(:) -1;
    U_wb_R_hi_wrap(end+1:end+length(U_wb_R_hi_bins(:))) = U_wb_R_hi_bins(:) +1;

    Dstroke_wb_hi_wrap(end+1:end+length(Dstroke_wb_hi_bins(:))) = Dstroke_wb_hi_bins(:) -1;
    Dstroke_wb_hi_wrap(end+1:end+length(Dstroke_wb_hi_bins(:))) = Dstroke_wb_hi_bins(:) +1;
    Dpitch_wb_hi_wrap(end+1:end+length(Dpitch_wb_hi_bins(:))) = Dpitch_wb_hi_bins(:) -1;
    Dpitch_wb_hi_wrap(end+1:end+length(Dpitch_wb_hi_bins(:))) = Dpitch_wb_hi_bins(:) +1;
    Ddev_wb_hi_wrap(end+1:end+length(Ddev_wb_hi_bins(:))) = Ddev_wb_hi_bins(:) -1;
    Ddev_wb_hi_wrap(end+1:end+length(Ddev_wb_hi_bins(:))) = Ddev_wb_hi_bins(:) +1;
    Daoa_wb_hi_wrap(end+1:end+length(Daoa_wb_hi_bins(:))) = Daoa_wb_hi_bins(:) -1;
    Daoa_wb_hi_wrap(end+1:end+length(Daoa_wb_hi_bins(:))) = Daoa_wb_hi_bins(:) +1;
    DU_wb_hi_wrap(end+1:end+length(DU_wb_hi_bins(:))) = DU_wb_hi_bins(:) -1;
    DU_wb_hi_wrap(end+1:end+length(DU_wb_hi_bins(:))) = DU_wb_hi_bins(:) +1;
    
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
    
    % L&R
    stroke_hi_func = csaps([t_hi_wrap;t_hi_wrap],[stroke_wb_L_hi_wrap;stroke_wb_R_hi_wrap]);
    pitch_hi_func = csaps([t_hi_wrap;t_hi_wrap],[pitch_wb_L_hi_wrap;pitch_wb_R_hi_wrap]);
    dev_hi_func = csaps([t_hi_wrap;t_hi_wrap],[dev_wb_L_hi_wrap;dev_wb_R_hi_wrap]);
    aoa_hi_func = csaps([t_hi_wrap;t_hi_wrap],[aoa_wb_L_hi_wrap;aoa_wb_R_hi_wrap]);
    U_hi_func = csaps([t_hi_wrap;t_hi_wrap],[U_wb_L_hi_wrap;U_wb_R_hi_wrap]);
end




if isempty(t_low_bins(:))==0
    t_low_wrap = t_low_bins(:);
    stroke_wb_L_low_wrap = stroke_wb_L_low_bins(:);
    stroke_wb_R_low_wrap = stroke_wb_R_low_bins(:);
    pitch_wb_L_low_wrap = pitch_wb_L_low_bins(:);
    pitch_wb_R_low_wrap = pitch_wb_R_low_bins(:);
    dev_wb_L_low_wrap = dev_wb_L_low_bins(:);
    dev_wb_R_low_wrap = dev_wb_R_low_bins(:);
    aoa_wb_L_low_wrap = aoa_wb_L_low_bins(:);
    aoa_wb_R_low_wrap = aoa_wb_R_low_bins(:);
    U_wb_L_low_wrap = U_wb_L_low_bins(:);
    U_wb_R_low_wrap = U_wb_R_low_bins(:);

    Dstroke_wb_low_wrap = Dstroke_wb_low_bins(:);
    Dpitch_wb_low_wrap = Dpitch_wb_low_bins(:);
    Ddev_wb_low_wrap = Ddev_wb_low_bins(:);
    Daoa_wb_low_wrap = Daoa_wb_low_bins(:);
    DU_wb_low_wrap = DU_wb_low_bins(:);

    t_low_wrap(end+1:end+length(t_low_bins(:))) = t_low_bins(:) -1;
    t_low_wrap(end+1:end+length(t_low_bins(:))) = t_low_bins(:) +1;
    
    stroke_wb_L_low_wrap(end+1:end+length(stroke_wb_L_low_bins(:))) = stroke_wb_L_low_bins(:) -1;
    stroke_wb_L_low_wrap(end+1:end+length(stroke_wb_L_low_bins(:))) = stroke_wb_L_low_bins(:) +1;
    pitch_wb_L_low_wrap(end+1:end+length(pitch_wb_L_low_bins(:))) = pitch_wb_L_low_bins(:) -1;
    pitch_wb_L_low_wrap(end+1:end+length(pitch_wb_L_low_bins(:))) = pitch_wb_L_low_bins(:) +1;
    dev_wb_L_low_wrap(end+1:end+length(dev_wb_L_low_bins(:))) = dev_wb_L_low_bins(:) -1;
    dev_wb_L_low_wrap(end+1:end+length(dev_wb_L_low_bins(:))) = dev_wb_L_low_bins(:) +1;
    aoa_wb_L_low_wrap(end+1:end+length(aoa_wb_L_low_bins(:))) = aoa_wb_L_low_bins(:) -1;
    aoa_wb_L_low_wrap(end+1:end+length(aoa_wb_L_low_bins(:))) = aoa_wb_L_low_bins(:) +1;
    U_wb_L_low_wrap(end+1:end+length(U_wb_L_low_bins(:))) = U_wb_L_low_bins(:) -1;
    U_wb_L_low_wrap(end+1:end+length(U_wb_L_low_bins(:))) = U_wb_L_low_bins(:) +1;
    
    stroke_wb_R_low_wrap(end+1:end+length(stroke_wb_R_low_bins(:))) = stroke_wb_R_low_bins(:) -1;
    stroke_wb_R_low_wrap(end+1:end+length(stroke_wb_R_low_bins(:))) = stroke_wb_R_low_bins(:) +1;
    pitch_wb_R_low_wrap(end+1:end+length(pitch_wb_R_low_bins(:))) = pitch_wb_R_low_bins(:) -1;
    pitch_wb_R_low_wrap(end+1:end+length(pitch_wb_R_low_bins(:))) = pitch_wb_R_low_bins(:) +1;
    dev_wb_R_low_wrap(end+1:end+length(dev_wb_R_low_bins(:))) = dev_wb_R_low_bins(:) -1;
    dev_wb_R_low_wrap(end+1:end+length(dev_wb_R_low_bins(:))) = dev_wb_R_low_bins(:) +1;
    aoa_wb_R_low_wrap(end+1:end+length(aoa_wb_R_low_bins(:))) = aoa_wb_R_low_bins(:) -1;
    aoa_wb_R_low_wrap(end+1:end+length(aoa_wb_R_low_bins(:))) = aoa_wb_R_low_bins(:) +1;
    U_wb_R_low_wrap(end+1:end+length(U_wb_R_low_bins(:))) = U_wb_R_low_bins(:) -1;
    U_wb_R_low_wrap(end+1:end+length(U_wb_R_low_bins(:))) = U_wb_R_low_bins(:) +1;

    Dstroke_wb_low_wrap(end+1:end+length(Dstroke_wb_low_bins(:))) = Dstroke_wb_low_bins(:) -1;
    Dstroke_wb_low_wrap(end+1:end+length(Dstroke_wb_low_bins(:))) = Dstroke_wb_low_bins(:) +1;
    Dpitch_wb_low_wrap(end+1:end+length(Dpitch_wb_low_bins(:))) = Dpitch_wb_low_bins(:) -1;
    Dpitch_wb_low_wrap(end+1:end+length(Dpitch_wb_low_bins(:))) = Dpitch_wb_low_bins(:) +1;
    Ddev_wb_low_wrap(end+1:end+length(Ddev_wb_low_bins(:))) = Ddev_wb_low_bins(:) -1;
    Ddev_wb_low_wrap(end+1:end+length(Ddev_wb_low_bins(:))) = Ddev_wb_low_bins(:) +1;
    Daoa_wb_low_wrap(end+1:end+length(Daoa_wb_low_bins(:))) = Daoa_wb_low_bins(:) -1;
    Daoa_wb_low_wrap(end+1:end+length(Daoa_wb_low_bins(:))) = Daoa_wb_low_bins(:) +1;
    DU_wb_low_wrap(end+1:end+length(DU_wb_low_bins(:))) = DU_wb_low_bins(:) -1;
    DU_wb_low_wrap(end+1:end+length(DU_wb_low_bins(:))) = DU_wb_low_bins(:) +1;
    
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
    
    % L&R
    stroke_low_func = csaps([t_low_wrap;t_low_wrap],[stroke_wb_L_low_wrap;stroke_wb_R_low_wrap]);
    pitch_low_func = csaps([t_low_wrap;t_low_wrap],[pitch_wb_L_low_wrap;pitch_wb_R_low_wrap]);
    dev_low_func = csaps([t_low_wrap;t_low_wrap],[dev_wb_L_low_wrap;dev_wb_R_low_wrap]);
    aoa_low_func = csaps([t_low_wrap;t_low_wrap],[aoa_wb_L_low_wrap;aoa_wb_R_low_wrap]);
    U_low_func = csaps([t_low_wrap;t_low_wrap],[U_wb_L_low_wrap;U_wb_R_low_wrap]);
end
