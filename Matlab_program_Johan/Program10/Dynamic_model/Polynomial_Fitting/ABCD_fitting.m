function [ ABCD_fit ] = ABCD_fitting( Maneuver_Matrix, dt )


    wb_t_squared = (36*dt)^2;


    ABCD_fit = {};

    [m_A_theta,b_A_theta,r_A_theta,sm_A_theta,sb_A_theta]   =   lsqfitgm(Maneuver_Matrix(:,14),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_A_phi,b_A_phi,r_A_phi,sm_A_phi,sb_A_phi]             =   lsqfitgm(Maneuver_Matrix(:,16),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_B_theta,b_B_theta,r_B_theta,sm_B_theta,sb_B_theta]   =   lsqfitgm(Maneuver_Matrix(:,17),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_B_phi,b_B_phi,r_B_phi,sm_B_phi,sb_B_phi]             =   lsqfitgm(Maneuver_Matrix(:,19),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_C_theta,b_C_theta,r_C_theta,sm_C_theta,sb_C_theta]   =   lsqfitgm(Maneuver_Matrix(:,20),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_C_phi,b_C_phi,r_C_phi,sm_C_phi,sb_C_phi]             =   lsqfitgm(Maneuver_Matrix(:,22),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_D_theta,b_D_theta,r_D_theta,sm_D_theta,sb_D_theta]   =   lsqfitgm(Maneuver_Matrix(:,23),wb_t_squared*Maneuver_Matrix(:,7));
    
    [m_D_phi,b_D_phi,r_D_phi,sm_D_phi,sb_D_phi]             =   lsqfitgm(Maneuver_Matrix(:,25),wb_t_squared*Maneuver_Matrix(:,7));

    
    
    ABCD_fit.m_A_theta = m_A_theta;
    ABCD_fit.b_A_theta = b_A_theta;
    ABCD_fit.m_A_phi = m_A_phi;
    ABCD_fit.b_A_phi = b_A_phi;
    
    ABCD_fit.m_B_theta = m_B_theta;
    ABCD_fit.b_B_theta = b_B_theta;
    ABCD_fit.m_B_phi = m_B_phi;
    ABCD_fit.b_B_phi = b_B_phi;
    
    ABCD_fit.m_C_theta = m_C_theta;
    ABCD_fit.b_C_theta = b_C_theta;
    ABCD_fit.m_C_phi = m_C_phi;
    ABCD_fit.b_C_phi = b_C_phi;
    
    ABCD_fit.m_D_theta = m_D_theta;
    ABCD_fit.b_D_theta = b_D_theta;
    ABCD_fit.m_D_phi = m_D_phi;
    ABCD_fit.b_D_phi = b_D_phi;

    
    figure()
    plot(Maneuver_Matrix(:,14),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min((Maneuver_Matrix(:,14))):((max(Maneuver_Matrix(:,14))-min(Maneuver_Matrix(:,14)))/1000):(max(Maneuver_Matrix(:,14)));
    plot(reg_int,b_A_theta+m_A_theta*reg_int,'r')
    title('point A vs strokeplane roll acc')
    xlabel('theta')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,16),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,16)):((max(Maneuver_Matrix(:,16))-min(Maneuver_Matrix(:,16)))/1000):max(Maneuver_Matrix(:,16));
    plot(reg_int,b_A_phi+m_A_phi*reg_int,'r')
    title('point A vs strokeplane roll acc')
    xlabel('phi')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,17),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,17)):((max(Maneuver_Matrix(:,17))-min(Maneuver_Matrix(:,17)))/1000):max(Maneuver_Matrix(:,17));
    plot(reg_int,b_B_theta+m_B_theta*reg_int,'r')
    title('point B vs strokeplane roll acc')
    xlabel('theta')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,19),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,19)):((max(Maneuver_Matrix(:,19))-min(Maneuver_Matrix(:,19)))/1000):max(Maneuver_Matrix(:,19));
    plot(reg_int,b_B_phi+m_B_phi*reg_int,'r')
    title('point B vs strokeplane roll acc')
    xlabel('phi')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,20),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,20)):((max(Maneuver_Matrix(:,20))-min(Maneuver_Matrix(:,20)))/1000):max(Maneuver_Matrix(:,20));
    plot(reg_int,b_C_theta+m_C_theta*reg_int,'r')
    title('point C vs strokeplane roll acc')
    xlabel('theta')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,22),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,22)):((max(Maneuver_Matrix(:,22))-min(Maneuver_Matrix(:,22)))/1000):max(Maneuver_Matrix(:,22));
    plot(reg_int,b_C_phi+m_C_phi*reg_int,'r')
    title('point C vs strokeplane roll acc')
    xlabel('phi')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,23),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,23)):((max(Maneuver_Matrix(:,23))-min(Maneuver_Matrix(:,23)))/1000):max(Maneuver_Matrix(:,23));
    plot(reg_int,b_D_theta+m_D_theta*reg_int,'r')
    title('point D vs strokeplane roll acc')
    xlabel('theta')
    ylabel('omega x dot')
    hold off
    
    figure()
    plot(Maneuver_Matrix(:,25),wb_t_squared*Maneuver_Matrix(:,7),'o')
    hold on
    reg_int = min(Maneuver_Matrix(:,25)):((max(Maneuver_Matrix(:,25))-min(Maneuver_Matrix(:,25)))/1000):max(Maneuver_Matrix(:,25));
    plot(reg_int,b_D_phi+m_D_phi*reg_int,'r')
    title('point D vs strokeplane roll acc')
    xlabel('phi')
    ylabel('omega x dot')
    hold off
    
end

