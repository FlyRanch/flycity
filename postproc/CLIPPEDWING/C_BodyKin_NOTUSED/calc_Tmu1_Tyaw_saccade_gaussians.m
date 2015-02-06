angle_pre_min = t_start;
angle_pre_max = t_stop;

angle_pre_ref  = t_plot(:);

%% Tyaw dynamics
angle_pre  = t_plot(:);
angle_post = Myaw_norm(:);

angle_pre  = angle_pre(angle_pre_ref>angle_pre_min & angle_pre_ref<angle_pre_max);
angle_post = angle_post(angle_pre_ref>angle_pre_min & angle_pre_ref<angle_pre_max);

angle_pre(isnan(angle_post))=[];
angle_post(isnan(angle_post))=[];

cf = createFit_gaussian_3rdOrder(angle_pre,angle_post,0);

figure
hold on
plot(angle_pre,angle_post,'.r')
plot(cf)

%% T_R dynamics
angle_pre  = t_plot(:);
angle_post = M_R_norm(:);

angle_pre  = angle_pre(angle_pre_ref>angle_pre_min & angle_pre_ref<angle_pre_max);
angle_post = angle_post(angle_pre_ref>angle_pre_min & angle_pre_ref<angle_pre_max);

angle_pre(isnan(angle_post))=[];
angle_post(isnan(angle_post))=[];

figure
plot(angle_pre,angle_post,'.')


