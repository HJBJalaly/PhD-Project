function CompareTime2(QQx,Coef_Opt,Time,Degree,Mode)

if(strcmp(Mode,'CostA'))
    rQ=Degree(1);
    SubplotNUM=1;
elseif(strcmp(Mode,'CostB'))
    nn=Degree(1);
    rQ=Degree(2) ;
    rU=Degree(3) ;
    SubplotNUM=2;
elseif(strcmp(Mode,'CostC'))
    nn=Degree(1);
    rQ=Degree(2) ;
    rU=Degree(3) ;
    SubplotNUM=2;    
end


%%  Recontruct Trajectories

%%%% Initial

Q1_X0=QQx(1,:);
    Q1_X0=InRangeShifter(Q1_X0);
Q2_X0=QQx(2,:);
    Q2_X0=InRangeShifter(Q2_X0);
Q3_X0=QQx(3,:);
    Q3_X0=InRangeShifter(Q3_X0);

%%%% Optimized
Coef_Opt_Q1=Coef_Opt(1:(rQ+1));
Coef_Opt_Q2=Coef_Opt((rQ+1)+1:2*(rQ+1));
Coef_Opt_Q3=Coef_Opt(2*(rQ+1)+1:3*(rQ+1));
Coef_Opt_D1Q1=Coef_Opt_Q1(1:end-1).*(rQ:-1:1);
Coef_Opt_D1Q2=Coef_Opt_Q2(1:end-1).*(rQ:-1:1);
Coef_Opt_D1Q3=Coef_Opt_Q3(1:end-1).*(rQ:-1:1);
Coef_Opt_D2Q1=Coef_Opt_D1Q1(1:end-1).*(rQ-1:-1:1);
Coef_Opt_D2Q2=Coef_Opt_D1Q2(1:end-1).*(rQ-1:-1:1);
Coef_Opt_D2Q3=Coef_Opt_D1Q3(1:end-1).*(rQ-1:-1:1);

Q1_Opt=polyval(Coef_Opt_Q1,Time);
    Q1_Opt=InRangeShifter(Q1_Opt);
Q2_Opt=polyval(Coef_Opt_Q2,Time);
    Q2_Opt=InRangeShifter(Q2_Opt);
Q3_Opt=polyval(Coef_Opt_Q3,Time);
    Q3_Opt=InRangeShifter(Q3_Opt);
D1Q1_Opt=polyval(Coef_Opt_D1Q1,Time);
D1Q2_Opt=polyval(Coef_Opt_D1Q2,Time);
D1Q3_Opt=polyval(Coef_Opt_D1Q3,Time);
D2Q1_Opt=polyval(Coef_Opt_D2Q1,Time);
D2Q2_Opt=polyval(Coef_Opt_D2Q2,Time);
D2Q3_Opt=polyval(Coef_Opt_D2Q3,Time);

QVal_Opt=[Q1_Opt;Q2_Opt;Q3_Opt];



figure('name',['Joints trajectory : '])
    subplot(3,1,1)
    plot(Time,rad2deg( Q1_X0),'linewidth',2,...
        'Color',[0.87058824300766 0.490196079015732 0],...
        'linestyle','--')
    hold on
    plot(Time,rad2deg(Q1_Opt),'linewidth',2,'color','b','linestyle','-')
    hold off
    set(gca,'FontWeight','bold','FontSize',12)
    title('Jonits Trajectory','FontWeight','bold','FontName','mwa_cmb10','FontSize',16);
    grid on
    ylabel('q_1 (deg)','fontsize',16,'FontName','mwa_cmb10');
    legend('Initial','Optimized','Orientation','horizontal')

    subplot(3,1,2)
    plot(Time,rad2deg(Q2_X0),'linewidth',2,...
        'Color',[0.87058824300766 0.490196079015732 0],...
        'linestyle','--')
    hold on
    plot(Time,rad2deg(Q2_Opt),'linewidth',2,'color','b','linestyle','-')
    hold off
    set(gca,'FontWeight','bold','FontSize',12)
    grid on
    ylabel('q_2 (deg)','fontsize',16,'FontName','mwa_cmb10');
    legend('Initial','Optimized','Orientation','horizontal')
    
    subplot(3,1,3)
    plot(Time,rad2deg(Q3_X0),'linewidth',2,...
        'Color',[0.87058824300766 0.490196079015732 0],...
        'linestyle','--')
    hold on
    plot(Time,rad2deg(Q3_Opt),'linewidth',2,'color','b','linestyle','-')
    hold off
    set(gca,'FontWeight','bold','FontSize',12)
    grid on
    xlabel('Time (s)','fontsize',16,'FontName','mwa_cmb10');
    ylabel('q_3 (deg)','fontsize',16,'FontName','mwa_cmb10');
    legend('Initial','Optimized','Orientation','horizontal')


%% EF

% Pos=[Xef;Yef];
% 
% RPos=L*[cos(Q1)+cos(Q1+Q2)+cos(Q1+Q2+Q3);
%         sin(Q1)+sin(Q1+Q2)+sin(Q1+Q2+Q3)];
% 
% RMS=sqrt(sum(sum((RPos-Pos).^2))*Tres/(Time(end)-Time(1)));
% 
% if(strcmp(ShowFlag,'Show'))
%     
%     figure('name',['WorkSapce : ',Name])
%         plot(Pos(1,:),Pos(2,:),'linewidth',2,'linestyle','-','color','b')
%         title('WorkSpace','FontWeight','bold','FontName','mwa_cmb10');
%         hold on
%         plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-','color','r')
%         hold off
%         xlabel('x (m)','fontsize',12,'FontName','mwa_cmb10')
%         ylabel('y (m)','fontsize',12,'FontName','mwa_cmb10')
% 
%     
%     AnimBot3DOF(Time,QVal',L);
% end
%% Trajectory

    

end