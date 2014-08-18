% figure('name','Torque vs Angle۲')
        subplot(3,1,1)
        plot(Q_Opt(1,1:49),Torque_Opt(1,1:49),'r','linewidth',2)
        title('Torque Angle Profile','FontWeight','bold')
        hold on
        plot(Q_Opt(1,49:100),Torque_Opt(1,49:100),'b','linewidth',2)
        
        plot(Q_Opt(1,1),Torque_Opt(1,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_1 (rad)','fontsize',12)
        ylabel('\tau_1','fontsize',14)
        hold off
        grid on

        subplot(3,1,2)
        plot(Q_Opt(2,1:49),Torque_Opt(2,1:49),'r','linewidth',2)
        hold on
        plot(Q_Opt(2,49:100),Torque_Opt(2,49:100),'b','linewidth',2)
        plot(Q_Opt(2,1),Torque_Opt(2,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_2 (rad)','fontsize',12)
        ylabel('\tau_2','fontsize',14)
        hold off
        grid on

        subplot(3,1,3)
        plot(Q_Opt(3,1:49),Torque_Opt(3,1:49),'r','linewidth',2)
        hold on
        plot(Q_Opt(3,49:100),Torque_Opt(3,49:100),'b','linewidth',2)
        
        plot(Q_Opt(3,1),Torque_Opt(3,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_3 (rad)','fontsize',12)
        ylabel('\tau_3','fontsize',14)
        hold off
        grid on
