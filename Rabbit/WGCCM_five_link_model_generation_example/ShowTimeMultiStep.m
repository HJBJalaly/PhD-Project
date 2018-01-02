function  ShowTimeMultiStep(Time,T_impact,Qq,DQq,MotionData,p_tib2,v_hip,mu,Ma)


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
    ylabel('q_1(rad)')
    xlabel('time(s)')
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
    ylabel('q_2(rad)')
    xlabel('time(s)')
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
    ylabel('q_3(rad)')
    xlabel('time(s)')
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
    ylabel('q_4(rad)')
    xlabel('time(s)')
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
    ylabel('q_5(rad)')
    xlabel('time(s)')
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
    subplot(3,2,1)
    plot(Time,DQq(:,1))
    ylabel('Dq_1(rad/s)')
    xlabel('time(s)')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,2)
    plot(Time,DQq(:,2))
    ylabel('Dq_2(rad/s)')
    xlabel('time(s)')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,3)
    plot(Time,DQq(:,3))
    ylabel('Dq_3(rad/s)')
    xlabel('time(s)')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,4)
    plot(Time,DQq(:,4))
    ylabel('Dq_4(rad/s)')
    xlabel('time(s)')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,5)
    plot(Time,DQq(:,5))
    ylabel('Dq_5(rad/s)')
    xlabel('time(s)')
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


figure(9)
    plot(Time,MotionData(5,:))% Ft
    hold all
    plot(Time,MotionData(6,:))%Fn
    plot(Time,MotionData(6,:)*mu,'linestyle','--','color','k')%Fn*mu
    plot(Time,-MotionData(6,:)*mu,'linestyle','--','color','k')%Fn*mu
    grid on
    legend('F_t','F_n','\mu*F_n')
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


figure(10)
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
    
figure(11)
    subplot(1,3,1)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(Qq(:,1),DQq(:,1),'linewidth',2)
    grid on
    hold on
    for i=2:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,1),DQq(index,1),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('\boldmath$q_1$ (deg)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    ylabel('\boldmath$\dot{q}_1$ (deg/s)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    
%     subplot(3,2,2)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(Qq(:,2),DQq(:,2),'linewidth',2)
%     grid on
%     hold on
%     for i=2:length(T_impact)
%         index=find(Time==T_impact(i),1);
%         plot(Qq(index,2),DQq(index,2),'color','r','linestyle','none','marker','*')
%     end
%     hold off
%     xlabel('\boldmath$q_2$ (deg)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
%     ylabel('\boldmath$\dot{q}_2$ (deg/s)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    
    
    subplot(1,3,2)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(Qq(:,3),DQq(:,3),'linewidth',2)
    grid on
    hold on
    for i=2:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,3),DQq(index,3),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('\boldmath$q_3$ (deg)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    ylabel('\boldmath$\dot{q}_3$ (deg/s)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    
%     subplot(3,2,4)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(Qq(:,4),DQq(:,4),'linewidth',2)
%     grid on
%     hold on
%     for i=2:length(T_impact)
%         index=find(Time==T_impact(i),1);
%         plot(Qq(index,4),DQq(index,4),'color','r','linestyle','none','marker','*')
%     end
%     hold off
%     xlabel('\boldmath$q_4$ (deg)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
%     ylabel('\boldmath$\dot{q}_4$ (deg/s)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    
    
    subplot(1,3,3)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(Qq(:,5),DQq(:,5),'linewidth',2)
    grid on
    hold on
    for i=2:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,5),DQq(index,5),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('\boldmath$q_5$ (deg)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    ylabel('\boldmath$\dot{q}_5$ (deg/s)','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    

figure(12)
    
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
    legend('Stance','Swing','orientation','horizontal')
    hold off
    subplot(3,2,2)
    legend('Swing','Stance','orientation','horizontal')
    hold off
    subplot(3,2,3)
    hold off
    subplot(3,2,4)
    hold off
    subplot(3,2,5)
    hold off
    subplot(3,2,6)
    hold off

figure(13)
    set(gcf,'position',[6 357 1266  195])
    i=length(T_impact)/2;
    indx1=find(Time==T_impact(2*(i)-2),1)+1;
    indx2=find(Time==T_impact(2*(i)-1),1);
    indx3=find(Time==T_impact(2*(i)-0),1);
    
    subplot(1,2,1)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    % for legend
    plot(rad2deg(Qq(indx1+100:indx1+100,1)),MotionData(1,indx1+100:indx1+100),'linewidth',2,'marker','o','markerfacecolor','b')
    hold on
    plot(rad2deg(Qq(indx2+1+100:indx2+1+100,2)),MotionData(2,indx2+1+100:indx2+1+100),'r','linewidth',2,'marker','^','markerfacecolor','r')
    legend('Stance','Swing','orientation','horizontal')
    %for marker
    plot(rad2deg(Qq(indx1+100:200:indx2,1)),MotionData(1,indx1+100:200:indx2),'linestyle','none','marker','o','markersize',8,'markerfacecolor','b')
    plot(rad2deg(Qq(indx2+1+100:200:indx3,2)),MotionData(2,indx2+1+100:200:indx3),'r','linestyle','none','marker','^','markersize',8,'markerfacecolor','r')
    % for  lines
    plot(rad2deg(Qq(indx1:indx2,1)),MotionData(1,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1:indx3,2)),MotionData(2,indx2+1:indx3),'r','linewidth',2)
    % for start point
    plot(rad2deg(Qq(indx1,1)),MotionData(1,indx1),'marker','*','markersize',8)
    plot(rad2deg(Qq(indx2+1,2)),MotionData(2,indx2+1),'r','marker','*','markersize',8)
    grid on
    hold off
    xlabel('\boldmath$q_1$ (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    ylabel('\boldmath$u_1$ (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    
%     subplot(3,2,2)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(Qq(indx1:indx2  ,2)),MotionData(2,indx1:indx2),'linewidth',2)
%     hold on
%     plot(rad2deg(Qq(indx2+1:indx3,1)),MotionData(1,indx2+1:indx3),'r','linewidth',2)
%     plot(rad2deg(Qq(indx1  ,2)),MotionData(2,indx1),'marker','*')
%     plot(rad2deg(Qq(indx2+1,1)),MotionData(1,indx2+1),'r','marker','*')
%     grid on
%     hold off
%     xlabel('q_2 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('u_2 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     legend('Swing','Stance','orientation','horizontal')
    
    subplot(1,2,2)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    % for legend
    plot(rad2deg(Qq(indx1+100:indx1+100,3)),MotionData(3,indx1+100:indx1+100),'linewidth',2,'marker','o','markerfacecolor','b')
    hold on
    plot(rad2deg(Qq(indx2+1+100:indx2+1+100,4)),MotionData(4,indx2+1+100:indx2+1+100),'r','linewidth',2,'marker','^','markerfacecolor','r')
    legend('Stance','Swing','orientation','horizontal')
    %for marker
    plot(rad2deg(Qq(indx1+100:200:indx2,3)),MotionData(3,indx1+100:200:indx2),'linestyle','none','marker','o','markersize',8,'markerfacecolor','b')
    plot(rad2deg(Qq(indx2+1+100:200:indx3,4)),MotionData(4,indx2+1+100:200:indx3),'r','linestyle','none','marker','^','markersize',8,'markerfacecolor','r')
    % for  lines
    plot(rad2deg(Qq(indx1:indx2  ,3)),MotionData(3,indx1:indx2),'linewidth',2)
    plot(rad2deg(Qq(indx2+1:indx3,4)),MotionData(4,indx2+1:indx3),'r','linewidth',2)
    % for start point
    plot(rad2deg(Qq(indx1 ,3)),MotionData(3,indx1),'marker','*','markersize',8)
    plot(rad2deg(Qq(indx2+1,4)),MotionData(4,indx2+1),'r','marker','*','markersize',8)
    grid on
    hold off
    xlabel('\boldmath$q_3$ (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10','interpret','latex');
    ylabel('\boldmath$u_3$ (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10','interpret','latex');

%     subplot(3,2,4)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(Qq(indx1  ,4)),MotionData(4,indx1),'marker','*')
%     hold on
%     plot(rad2deg(Qq(indx1:indx2  ,4)),MotionData(4,indx1:indx2),'linewidth',2)
%     plot(rad2deg(Qq(indx2+1,3)),MotionData(3,indx2+1),'r','marker','*')
%     plot(rad2deg(Qq(indx2+1:indx3,3)),MotionData(3,indx2+1:indx3),'r','linewidth',2)
%     grid on
%     hold off
%     xlabel('q_4 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('u_4 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
% 
%     subplot(3,2,5)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(Qq(indx1 ,5)),MotionData(5,indx1),'marker','*')
%     hold on
%     plot(rad2deg(Qq(indx1:indx2  ,5)),MotionData(5,indx1:indx2),'linewidth',2)
%     plot(rad2deg(Qq(indx2+1,5)),MotionData(5,indx2+1),'r','marker','*')
%     plot(rad2deg(Qq(indx2+1:indx3,5)),MotionData(5,indx2+1:indx3),'r','linewidth',2)
%     grid on
%     hold off
%     xlabel('q_5 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('u_5 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
   
figure(14)    
    i=length(T_impact)/2;
    indx1=find(Time==T_impact(2*(i)-2),1)+1;
    indx2=find(Time==T_impact(2*(i)-1),1);
    indx3=find(Time==T_impact(2*(i)-0),1);
    
    subplot(1,2,1)
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    plot(rad2deg(DQq(indx1:indx2,1)),MotionData(1,indx1:indx2),'linewidth',2)
    hold on
    plot(rad2deg(DQq(indx2+1:indx3,2)),MotionData(2,indx2+1:indx3),'r','linewidth',2)
    plot(rad2deg(DQq(indx1,1)),MotionData(1,indx1),'marker','*')
    plot(rad2deg(DQq(indx2+1,2)),MotionData(2,indx2+1),'r','marker','*')
    grid on
    hold off
    xlabel('Dq_1 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    legend('Stance','Swing','orientation','horizontal')

%     subplot(3,2,2)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(DQq(indx1:indx2  ,2)),MotionData(2,indx1:indx2),'linewidth',2)
%     hold on
%     plot(rad2deg(DQq(indx2+1:indx3,1)),MotionData(1,indx2+1:indx3),'r','linewidth',2)
%     plot(rad2deg(DQq(indx1  ,2)),MotionData(2,indx1),'marker','*')
%     plot(rad2deg(DQq(indx2+1,1)),MotionData(1,indx2+1),'r','marker','*')
%     grid on
%     hold off
%     xlabel('Dq_2 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('u_2 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     legend('Swing','Stance','orientation','horizontal')
    
    subplot(1,2,2)
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

%     subplot(3,2,4)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(DQq(indx1  ,4)),MotionData(4,indx1),'marker','*')
%     hold on
%     plot(rad2deg(DQq(indx1:indx2  ,4)),MotionData(4,indx1:indx2),'linewidth',2)
%     plot(rad2deg(DQq(indx2+1,3)),MotionData(3,indx2+1),'r','marker','*')
%     plot(rad2deg(DQq(indx2+1:indx3,3)),MotionData(3,indx2+1:indx3),'r','linewidth',2)
%     grid on
%     hold off
%     xlabel('Dq_4 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('u_4 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

%     subplot(3,2,5)
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(DQq(indx1 ,5)),MotionData(5,indx1),'marker','*')
%     hold on
%     plot(rad2deg(DQq(indx1:indx2  ,5)),MotionData(5,indx1:indx2),'linewidth',2)
%     plot(rad2deg(DQq(indx2+1,5)),MotionData(5,indx2+1),'r','marker','*')
%     plot(rad2deg(DQq(indx2+1:indx3,5)),MotionData(5,indx2+1:indx3),'r','linewidth',2)
%     grid on
%     hold off
%     xlabel('Dq_5 (deg/s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('u_5 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    
figure(15)
    subplot(2,1,1)
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time(indx1:indx2),MotionData(1,indx1:indx2),'linewidth',2)
        hold on
        plot(Time(indx2+1:indx3),MotionData(2,indx2+1:indx3),'r','linewidth',2)
        hold off
        grid on
        xlabel('Time (s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('u_1 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        legend('Stance','Swing','orientation','horizontal')
    subplot(2,1,2)
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time(indx1:indx2),MotionData(3,indx1:indx2),'linewidth',2)
        hold on
        plot(Time(indx2+1:indx3),MotionData(4,indx2+1:indx3),'r','linewidth',2)
        hold off
        grid on
        xlabel('Time (s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('u_3 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

figure(16)
    subplot(2,1,1)
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time(indx1:indx2),MotionData(1,indx1:indx2)'.*DQq(indx1:indx2  ,1),'linewidth',2)
        hold on
        plot(Time(indx2+1:indx3),MotionData(2,indx2+1:indx3)'.*DQq(indx2+1:indx3  ,2),'r','linewidth',2)
        hold off
        grid on
        xlabel('Time (s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('power_h_i_p (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        legend('Stance','Swing','orientation','horizontal')
    subplot(2,1,2)
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time(indx1:indx2),MotionData(3,indx1:indx2)'.*DQq(indx1:indx2  ,3),'r','linewidth',2)
        hold on
        plot(Time(indx2+1:indx3),MotionData(4,indx2+1:indx3)'.*DQq(indx2+1:indx3  ,4),'linewidth',2)
        grid on
        hold off
        xlabel('Time (s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('power_k_n_e_e (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');


figure(17)
    subplot(2,1,1)
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time(indx1:indx2),rad2deg(Qq(indx1:indx2  ,1)),'linewidth',2)
        hold on
        plot(Time(indx2+1:indx3),rad2deg(Qq(indx2+1:indx3  ,2)),'r','linewidth',2)
        grid on
        hold off
        xlabel('Time (s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('q_1 (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        legend('Stance','Swing','orientation','horizontal')
    subplot(2,1,2)
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time(indx1:indx2),rad2deg(Qq(indx1:indx2  ,3)),'linewidth',2)
        hold on
        plot(Time(indx2+1:indx3),rad2deg(Qq(indx2+1:indx3  ,4)),'r','linewidth',2)
        grid on
        hold off
        xlabel('Time (s)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('q_3 (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');


    