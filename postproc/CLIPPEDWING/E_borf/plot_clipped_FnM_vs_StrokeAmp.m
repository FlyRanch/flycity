            
            % F-Amp
%             figure
            subplot(1,2,1)
            hold off
            plot(Amp_ratio(n_now),Fx_norm(n_now),'ok','markerfacecolor','b')
            hold on
            plot(Amp_ratio(n_now),Fy_norm_MinSteady(n_now),'ok','markerfacecolor','r')
            plot(Amp_ratio(n_now),Fz_norm(n_now),'ok','markerfacecolor','g')

            % linear fits
            plot([min(Amp_ratio(n_now)) max(Amp_ratio(n_now))],polyval(Fx_Amp_fit(p,:),[min(Amp_ratio(n_now)) max(Amp_ratio(n_now))]),'k')
            plot([min(Amp_ratio(n_now)) max(Amp_ratio(n_now))],polyval(Fy_MinSteady_Amp_fit(p,:),[min(Amp_ratio(n_now)) max(Amp_ratio(n_now))]),'k')
            plot([min(Amp_ratio(n_now)) max(Amp_ratio(n_now))],polyval(Fz_Amp_fit(p,:),[min(Amp_ratio(n_now)) max(Amp_ratio(n_now))]),'k')

            legend('x','y','z','location','E')
            xlabel('Amp ratio')
            ylabel('F/mg')
%             axis([1 1.25 floor(10*min([Fx_norm(n_now);Fy_norm_MinSteady(n_now);Fz_norm(n_now)]))/10 ceil(10*max([Fx_norm(n_now);Fy_norm_MinSteady(n_now);Fz_norm(n_now)]))/10])
            axis([1 1.25 -1.5 .25])
            set(gca,'xtick',0.75:.25:1.25)
            set(gca,'ytick',-1.5:.25:1.5)

            % M-Amp
            subplot(1,2,2)
            hold off
            plot(Amp_ratio(n_now),Mx_norm_MinSteady(n_now),'ok','markerfacecolor','b')
            hold on
            plot(Amp_ratio(n_now),My_norm_CoM(n_now),'ok','markerfacecolor','r')
            plot(Amp_ratio(n_now),Mz_norm_MinSteady(n_now),'ok','markerfacecolor','g')

            plot([min(Amp_ratio(n_now)) max(Amp_ratio(n_now))],polyval(Mx_MinSteady_Amp_fit(p,:),[min(Amp_ratio(n_now)) max(Amp_ratio(n_now))]),'k')
            plot([min(Amp_ratio(n_now)):.01:max(Amp_ratio(n_now))],polyval(My_CoM_Amp_fit2(p,:),[min(Amp_ratio(n_now)):.01:max(Amp_ratio(n_now))]),'k')
            plot([min(Amp_ratio(n_now)) max(Amp_ratio(n_now))],polyval(Mz_MinSteady_Amp_fit(p,:),[min(Amp_ratio(n_now)) max(Amp_ratio(n_now))]),'k')

%             legend('Mx','My@CoM','Mz','location','SW')
            xlabel('Amp ratio')
            ylabel('T/mgl')
%             axis([1 1.25 floor(10*min([Mx_norm_MinSteady(n_now);My_norm_CoM(n_now);Mz_norm_MinSteady(n_now)]))/10 ceil(10*max([Mx_norm_MinSteady(n_now);My_norm_CoM(n_now);Mz_norm_MinSteady(n_now)]))/10])
            axis([1 1.25 -.3 .3])
            set(gca,'xtick',0.75:.25:1.25)
            set(gca,'ytick',-1.5:.1:1.5)
