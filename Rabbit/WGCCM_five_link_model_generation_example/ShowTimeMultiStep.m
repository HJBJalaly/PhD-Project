function  ShowTimeMultiStep(Time,T_impact,Qq,DQq,MotionData,p_tib2,v_hip)


figure(6)
    plot(Time,MotionData(1:4,:)','linewidth',2)
    legend('u_1','u_2','u_3','u_4')
    grid on
    xlabel('Time (s)')
    ylabel('Torque (N.m)')
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','k','linestyle','-.','linewidth',1)   
    end
    hold off
    
    

figure(7)
    subplot(3,2,1)
    plot(Time,Qq(:,1))
    ylabel('q_1')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,2)
    plot(Time,Qq(:,2))
    ylabel('q_2')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,3)
    plot(Time,Qq(:,3))
    ylabel('q_3')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,4)
    plot(Time,Qq(:,4))
    ylabel('q_4')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,5)
    plot(Time,Qq(:,5))
    ylabel('q_5')
    grid on
    hold on
    h=get(gca,'YLIM')';
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
        drawnow
    end
    drawnow
    hold off
    
    subplot(3,2,6)
    plot(0,0)
    

figure(8)
    plot(Time,MotionData(5,:))
    hold all
    plot(Time,MotionData(6,:))
    grid on
    legend('F_t','F_n')
    hold off
    xlabel('Time (s)')
    ylabel('GRF (N)')
        h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.','linewidth',2)   
    end
    hold off


figure(9)
    subplot(2,1,1)
    plot(Time,p_tib2(1,:))
    hold all
    plot(Time,p_tib2(2,:))
    grid on
    legend('p^x_t_2','p^y_t_2')
    hold off
    xlabel('Time (s)')
    ylabel('p_t_i_b_2 (m)')
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.','linewidth',2)   
    end
    hold off
    
    subplot(2,1,2)
    plot(Time,v_hip(1,:))
    hold all
    plot(Time,v_hip(2,:))
    grid on
    legend('v^x_h','v^y_h')
    hold off
    xlabel('Time (s)')
    ylabel('v_h_i_p (m/s)')
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.','linewidth',2)   
    end
    hold off
    
figure(10)
    subplot(3,2,1)
    plot(Qq(:,1),DQq(:,1))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,1),DQq(index,1),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_1')
    ylabel('dq_1')
    
    subplot(3,2,2)
    plot(Qq(:,2),DQq(:,2))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,2),DQq(index,2),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_2')
    ylabel('dq_2')
    
    
    subplot(3,2,3)
    plot(Qq(:,3),DQq(:,3))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,3),DQq(index,3),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_3')
    ylabel('dq_3')
    
    subplot(3,2,4)
    plot(Qq(:,4),DQq(:,4))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,4),DQq(index,4),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_4')
    ylabel('dq_4')
    
    subplot(3,2,5)
    plot(Qq(:,5),DQq(:,5))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,5),DQq(index,5),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_5')
    ylabel('dq_5')

figure(11)
    
    indx1=1;
    for i=1:length(T_impact)/2
        indx2=find(Time==T_impact(2*(i)-1),1);
        indx3=find(Time==T_impact(2*(i)-0),1);
    
        subplot(3,2,1)
        plot(Qq(indx1:indx2,1),MotionData(1,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,2),MotionData(2,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_1 without relabling (rad)')
        ylabel('u_1 without relabling (N.m)')

        subplot(3,2,2)
        plot(Qq(indx1:indx2  ,2),MotionData(2,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,1),MotionData(1,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_2 without relabling (rad)')
        ylabel('u_2 without relabling (N.m)')


        subplot(3,2,3)
        plot(Qq(indx1:indx2  ,3),MotionData(3,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,4),MotionData(4,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_3 without relabling (rad)')
        ylabel('u_3 without relabling (N.m)')

        subplot(3,2,4)
        plot(Qq(indx1:indx2  ,4),MotionData(4,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,3),MotionData(3,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_4 without relabling (rad)')
        ylabel('u_4 without relabling (N.m)')

        subplot(3,2,5)
        plot(Qq(indx1:indx2  ,5),MotionData(5,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,5),MotionData(5,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_5 without relabling (rad)')
        ylabel('u_5 without relabling (N.m)')

       indx1=indx2+1;
    end
    subplot(3,2,1)
    hold off
    subplot(3,2,2)
    hold off
    subplot(3,2,3)
    hold off
    subplot(3,2,4)
    hold off
    subplot(3,2,5)
    hold off
    subplot(3,2,6)
    hold off

figure(12)    
    i=length(T_impact)/2;
    indx1=find(Time==T_impact(2*(i)-2),1)+1;
    indx2=find(Time==T_impact(2*(i)-1),1);
    indx3=find(Time==T_impact(2*(i)-0),1);
    
    subplot(3,2,1)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(Qq(indx1,1)),MotionData(1,indx1),'marker','*')
    hold on
    plot(rad2deg(Qq(indx1:indx2,1)),MotionData(1,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1,2)),MotionData(2,indx2+1),'r','marker','*')
    plot(rad2deg(Qq(indx2+1:indx3,2)),MotionData(2,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    subplot(3,2,2)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(Qq(indx1  ,2)),MotionData(2,indx1),'marker','*')
    hold on
    plot(rad2deg(Qq(indx1:indx2  ,2)),MotionData(2,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1,1)),MotionData(1,indx2+1),'r','marker','*')
    plot(rad2deg(Qq(indx2+1:indx3,1)),MotionData(1,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');


    subplot(3,2,3)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(Qq(indx1 ,3)),MotionData(3,indx1),'marker','*')
    hold on
    plot(rad2deg(Qq(indx1:indx2  ,3)),MotionData(3,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1,4)),MotionData(4,indx2+1),'r','marker','*')
    plot(rad2deg(Qq(indx2+1:indx3,4)),MotionData(4,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('q_3 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    subplot(3,2,4)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(Qq(indx1  ,4)),MotionData(4,indx1),'marker','*')
    hold on
    plot(rad2deg(Qq(indx1:indx2  ,4)),MotionData(4,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1,3)),MotionData(3,indx2+1),'r','marker','*')
    plot(rad2deg(Qq(indx2+1:indx3,3)),MotionData(3,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('q_4 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_4 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    subplot(3,2,5)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(Qq(indx1 ,5)),MotionData(5,indx1),'marker','*')
    hold on
    plot(rad2deg(Qq(indx1:indx2  ,5)),MotionData(5,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1,5)),MotionData(5,indx2+1),'r','marker','*')
    plot(rad2deg(Qq(indx2+1:indx3,5)),MotionData(5,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('q_5 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_5 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

   
figure(13)    
    i=length(T_impact)/2;
    indx1=find(Time==T_impact(2*(i)-2),1)+1;
    indx2=find(Time==T_impact(2*(i)-1),1);
    indx3=find(Time==T_impact(2*(i)-0),1);
    
    subplot(3,2,1)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(DQq(indx1,1)),MotionData(1,indx1),'marker','*')
    hold on
    plot(rad2deg(DQq(indx1:indx2,1)),MotionData(1,indx1:indx2),'linewidth',2)
    plot(rad2deg(DQq(indx2+1,2)),MotionData(2,indx2+1),'r','marker','*')
    plot(rad2deg(DQq(indx2+1:indx3,2)),MotionData(2,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('Dq_1 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    subplot(3,2,2)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(DQq(indx1  ,2)),MotionData(2,indx1),'marker','*')
    hold on
    plot(rad2deg(DQq(indx1:indx2  ,2)),MotionData(2,indx1:indx2),'linewidth',2)
    plot(rad2deg(DQq(indx2+1,1)),MotionData(1,indx2+1),'r','marker','*')
    plot(rad2deg(DQq(indx2+1:indx3,1)),MotionData(1,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('Dq_2 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');


    subplot(3,2,3)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(DQq(indx1 ,3)),MotionData(3,indx1),'marker','*')
    hold on
    plot(rad2deg(DQq(indx1:indx2  ,3)),MotionData(3,indx1:indx2),'linewidth',2)
    plot(rad2deg(DQq(indx2+1,4)),MotionData(4,indx2+1),'r','marker','*')
    plot(rad2deg(DQq(indx2+1:indx3,4)),MotionData(4,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('Dq_3 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    subplot(3,2,4)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(DQq(indx1  ,4)),MotionData(4,indx1),'marker','*')
    hold on
    plot(rad2deg(DQq(indx1:indx2  ,4)),MotionData(4,indx1:indx2),'linewidth',2)
    plot(rad2deg(DQq(indx2+1,3)),MotionData(3,indx2+1),'r','marker','*')
    plot(rad2deg(DQq(indx2+1:indx3,3)),MotionData(3,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('Dq_4 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_4 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    subplot(3,2,5)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(DQq(indx1 ,5)),MotionData(5,indx1),'marker','*')
    hold on
    plot(rad2deg(DQq(indx1:indx2  ,5)),MotionData(5,indx1:indx2),'linewidth',2)
    plot(rad2deg(DQq(indx2+1,5)),MotionData(5,indx2+1),'r','marker','*')
    plot(rad2deg(DQq(indx2+1:indx3,5)),MotionData(5,indx2+1:indx3),'r','linewidth',2)
    grid on
    hold off
    xlabel('Dq_5 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_5 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
