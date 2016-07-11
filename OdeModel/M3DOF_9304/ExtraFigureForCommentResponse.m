%%
load testCc_eli_2_rb0
% load testCc_Lin_1_WithOptimizedTrajectory3_ForPaper_WindowsOnly_Uni_rQ_20
% figure(1)
% clf
figure(2)
clf
ActuationCost=[];
ComplexCost=[];
for i=-2:2:2
    
    Landa(:)=2e-3*(10)^i
    Landa(1)=0;
    

    [TorqueDesire_Opt,TorqueActive_Opt,TorquePassiveOptimal,TorqueBicepsOptimal,Q_Opt,D1Q_Opt,D2Q_Opt,BetaOptimal_Opt,ThetaOptimal_Opt,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,IntAbsUdqDesire_Opt,CostSlopeD1Q_Opt,CostSlopeD2Q_Opt,CostParam_Opt,RMSError_Opt]=...
                           ShowTime(x      ,Time,Tres,Degree,Weight,Landa,SampleRate,[],[],[],XEF,YEF,m,L,g,[],'DntShow','1Cycle','CostCc','Optimized');

    ActuationCost(end+1)=(IntU2_Opt)*Tres;
    ComplexCost(end+1)=CostParam_Opt(1)*Tres;


    Q1=Q_Opt(1,:);
    Q2=Q_Opt(2,:);
    Q3=Q_Opt(3,:);
    
%     figure(1)
%     subplot(5,3,(i+2)*3+1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
%     plot(rad2deg(Q1),TorquePassiveOptimal(1,:),'linewidth',3,'color','g')
%     xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     grid on
%     set(gca,'YMinorGrid','on')
%     hold on
%     plot(rad2deg( Q1),TorqueDesire(1,:),'linewidth',2,'linestyle','-.')
%     hold off
%     legend('u_m_1','u_r_1','Orientation','horizontal')
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     %
%     subplot(5,3,(i+2)*3+2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(Q2),TorquePassiveOptimal(2,:),'linewidth',3,'color','g')
%     xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     grid on
%     set(gca,'YMinorGrid','on')
%     hold on
%     plot(rad2deg( Q2),TorqueDesire(2,:),'linewidth',2,'linestyle','-.')
%     hold off
%     legend('u_m_2','u_r_2','Orientation','horizontal')
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     %
%     subplot(5,3,(i+2)*3+3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     plot(rad2deg(Q3),TorquePassiveOptimal(3,:),'linewidth',3,'color','g')
%     xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     grid on
%     set(gca,'YMinorGrid','on')
%     hold on
%     plot(rad2deg(Q3),TorqueDesire(3,:),'linewidth',2,'linestyle','-.')
%     hold off
%     legend('u_m_3','u_r_3','Orientation','horizontal')
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     
    figure(2)
    subplot(3,1,1)
    plot(rad2deg(Q1),TorquePassiveOptimal(1,:),'linewidth',2)
    hold all
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_m_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    set(gca,'YMinorGrid','on')
    
    subplot(3,1,2)
    plot(rad2deg(Q2),TorquePassiveOptimal(2,:),'linewidth',3)
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_m_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    set(gca,'YMinorGrid','on')
    hold all
    subplot(3,1,3)
    plot(rad2deg(Q3),TorquePassiveOptimal(3,:),'linewidth',3)
    hold all
    xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_m_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    set(gca,'YMinorGrid','on')    
    
end
cd
%% 
load testCc_eli_3(ForPaper_JustOnWindows)
% Landa(2)=1e-3;
% Landa(4)=1e-3;
ActuationCost=[];
ComplexCost=[];

for i=[-4,0,4]
    
    Landa(1)=2e-4*(10)^i;
    Landa(3)=2e-4*(10)^i;

    [TorqueDesire,TorqueActive,TorquePassiveOptimal,TorqueBicepsOptimal,Q_Opt,D1Q_Opt,D2Q_Opt,BetaOptimal_Opt,ThetaOptimal_Opt,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,IntAbsUdqDesire_Opt,CostSlopeD1Q_Opt,CostSlopeD2Q_Opt,CostParam_Opt,RMSError_Opt]=...
                           ShowTime(x      ,Time,Tres,Degree,Weight,Landa,SampleRate,[],[],[],XEF,YEF,m,L,g,[],'DntShow','1Cycle','CostCc','Optimized');

    
    ActuationCost(end+1)=(IntU2_Opt)*Tres;
    ComplexCost(end+1)=CostParam_Opt(1)*Tres;
                   
    figure('name',['gamma=',num2str(Landa(1))])
    subplot(3,1,1)
%     plot(Time,TorqueDesire(1,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimal(1,:),'linewidth',2,'linestyle','-','color','g')
    plot(Time,TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-','Color',[0.75 0 0.75])
%     plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_m_1','u_b_1');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    %
    subplot(3,1,2)
%     plot(Time,TorqueDesire(2,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimal(2,:),'linewidth',2,'linestyle','-','color','g')
    plot(Time,TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-','Color',[0.75 0 0.75])
    plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','Color',[0.68 0.45 0])
%     plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_m_2','u_b_1','u_b_2');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    %
    subplot(3,1,3)
%     plot(Time,TorqueDesire(3,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimal(3,:),'linewidth',2,'linestyle','-','color','g')
    plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','Color',[0.68 0.45 0])
%     plot(Time,TorqueActive(3,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_m_3','u_b_2');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    xlabel('time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');

end