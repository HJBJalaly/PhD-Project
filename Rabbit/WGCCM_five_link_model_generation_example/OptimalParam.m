function  [CostActuation,CostD2Q,CostParaReg,BetaOptimal,ThetaOptimal,EtaOptimal,TorqueMonoSpring1,TorqueMonoSpring3,TorqueDamper1,TorqueDamper3,TorqueBiSpring13]=...
                    OptimalParam(Time,Q1,Q3,Qhat13,DtQ1,DtQ3,Ur1,Ur3,rM,rB,rD,Landa,Gamma,Weight,SampleRate,SwitchIndex,Show,estimate)

W1=Weight(1);
W3=Weight(2);
                
QQ_1=[];
QQ_3=[];
QQc_1=[];
QQc_3=[];
for i=1:SampleRate:length(Q1)
    QQ_1(end+1,:) = Q1(i).^(rM:-1:0)';
    QQ_3(end+1,:) = Q3(i).^(rM:-1:0)';
end
for i=1:length(Q1)
    QQc_1(end+1,:) = Q1(i).^(rM:-1:0)';
    QQc_3(end+1,:) = Q3(i).^(rM:-1:0)';
end


if(rD>=0)
    DD_1=[];
    DD_3=[];
    DDc_1=[];
    DDc_3=[];
    for i=1:SampleRate:length(DtQ1)
        DD_1(end+1,:) = DtQ1(i).^(2*rD+1:-2:1)';
        DD_3(end+1,:) = DtQ3(i).^(2*rD+1:-2:1)';
    end
    for i=1:length(DtQ1)
        DDc_1(end+1,:) = DtQ1(i).^(2*rD+1:-2:1)';
        DDc_3(end+1,:) = DtQ3(i).^(2*rD+1:-2:1)';
    end
    
else
    DDc_1=zeros(size(DtQ1));
    DDc_3=zeros(size(DtQ3));
end

if(rB>0)
    QH_13=[];
    QHc_13=[];
    for i=1:SampleRate:length(Qhat13)
        QH_13(end+1,:) = Qhat13(i).^(rB:-1:0)';
    end
    for i=1:length(Qhat13)
        QHc_13(end+1,:) = Qhat13(i).^(rB:-1:0)';
    end
    
else
    QHc_13=zeros(size(Qhat13));
end


A=[];
B=[];
D=[];
Y=[];

rMp=rM+1;
rBp=rB+1;
rDp=rD+1;

if(rM>0 && rD==-1 && rB==0) % mono

    i=1;
    A((i-1)*rMp+1 :(i-0)*rMp , (i-1)*rMp+1:(i-0)*rMp )=inv((1+Gamma)*W1*(QQ_1'*QQ_1)+Landa*eye(rMp));
    Y((i-1)*rMp+1 :(i-0)*rMp , 1 )=W1*(QQ_1'*Ur1(1:SampleRate:end));

    i=2;
    A((i-1)*rMp+1 :(i-0)*rMp , (i-1)*rMp+1:(i-0)*rMp )=inv((1+Gamma)*W3*(QQ_3'*QQ_3)+Landa*eye(rMp));
    Y((i-1)*rMp+1 :(i-0)*rMp , 1 )=W3*(QQ_3'*Ur3(1:SampleRate:end));
    
    Ai=A;
    OptimalParametrs=Ai*Y;

    BetaOptimal1  =OptimalParametrs(0*rMp+1 : 1*rMp );
    BetaOptimal3  =OptimalParametrs(1*rMp+1 : 2*rMp );
    ThetaOptimal13=0;

    EtaOptimal1   =0;
    EtaOptimal3   =0;
    
    
elseif(rM>0 && rD>-1 && rB==0) % mono+damper
    i=1;
    A((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )=inv([ (1+Gamma)*W1*(QQ_1'*QQ_1) , W1*(QQ_1'*DD_1);  W1*(DD_1'*QQ_1) ,(1+Gamma)*W1*(DD_1'*DD_1)]+Landa*eye(rMp+rDp));
    Amain((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )=([ (1+Gamma)*W1*(QQ_1'*QQ_1) , W1*(QQ_1'*DD_1);  W1*(DD_1'*QQ_1) ,(1+Gamma)*W1*(DD_1'*DD_1)]+Landa*eye(rMp+rDp));
    Y((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , 1 )=W1*([QQ_1';DD_1']*Ur1(1:SampleRate:end));

    i=2;
    A((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )=inv([ (1+Gamma)*W3*(QQ_3'*QQ_3) , W3*(QQ_3'*DD_3);  W3*(DD_3'*QQ_3) , (1+Gamma)*W3*(DD_3'*DD_3)]+Landa*eye(rMp+rDp));
    Amain((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )=([ (1+Gamma)*W3*(QQ_3'*QQ_3) , W3*(QQ_3'*DD_3);  W3*(DD_3'*QQ_3) , (1+Gamma)*W3*(DD_3'*DD_3)]+Landa*eye(rMp+rDp));
    Y((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , 1 )=W3*([QQ_3';DD_3']*Ur3(1:SampleRate:end));

    Ai=A;
    OptimalParametrs=Ai*Y;

    BetaOptimal1  =OptimalParametrs(0*rMp+0*rDp+1 : 1*rMp+0*rDp    );
    BetaOptimal3  =OptimalParametrs(1*rMp+1*rDp+1 : 2*rMp+1*rDp    );
    ThetaOptimal13=0;

    EtaOptimal1   =OptimalParametrs(1*rMp+0*rDp+1 : 1*rMp+1*rDp    );
    EtaOptimal1 =diag(1-sign(EtaOptimal1))/2*EtaOptimal1;
    EtaOptimal3   =OptimalParametrs(2*rMp+1*rDp+1 : 2*rMp+2*rDp    );
    EtaOptimal3 =diag(1-sign(EtaOptimal3))/2*EtaOptimal3;
    
    if(~estimate)
        H=Amain;    
        options = optimoptions('quadprog',...
            'Algorithm','interior-point-convex','Display','off','MaxIter',100,...
            'TolFun',.00001);
        Xstar=quadprog(H,-Y,[],[],[],[],-inf*ones(2*rMp+2*rDp,1),[inf*ones(rMp,1);zeros(rDp,1);inf*ones(rMp,1);zeros(rDp,1)],OptimalParametrs,options);
    
        BetaOptimal1  =Xstar(0*rMp+0*rDp+1 : 1*rMp+0*rDp    );
        BetaOptimal3  =Xstar(1*rMp+1*rDp+1 : 2*rMp+1*rDp    );
        ThetaOptimal13=0;

        EtaOptimal1   =Xstar(1*rMp+0*rDp+1 : 1*rMp+1*rDp    );
        EtaOptimal3   =Xstar(2*rMp+1*rDp+1 : 2*rMp+2*rDp    );
    end

    

elseif(rM>0 && rD==-1 && rB>0) % mono+bi
    i=1;
    A((i-1)*rMp+1 :(i-0)*rMp , (i-1)*rMp+1:(i-0)*rMp )=inv((1+Gamma)*W1*(QQ_1'*QQ_1) +Landa*eye(rMp));
    B((i-1)*rMp+1 :(i-0)*rMp , 1:rBp )=W1*(QQ_1'*QH_13);
    Y((i-1)*rMp+1 :(i-0)*rMp , 1 )=W1*(QQ_1'*Ur1(1:SampleRate:end));

    i=2;
    A((i-1)*rMp+1 :(i-0)*rMp , (i-1)*rMp+1:(i-0)*rMp )=inv((1+Gamma)*W1*(QQ_3'*QQ_3) +Landa*eye(rMp));
    B((i-1)*rMp+1 :(i-0)*rMp , 1:rBp )=W1*(QQ_3'*QH_13);
    Y((i-1)*rMp+1 :(i-0)*rMp , 1 )=W1*(QQ_3'*Ur3(1:SampleRate:end));

    D=(W1+W3)*(1+Gamma)*(QH_13'*QH_13)+Landa*eye(rBp);
    Y((2)*rMp+1 :(2)*rMp+rBp , 1 )=W1*(QH_13'*Ur1(1:SampleRate:end))+W3*(QH_13'*Ur3(1:SampleRate:end));
    
    Ai=A;
    Delta=inv(D-B'*Ai*B);
    OptimalParametrs=[ Ai+Ai*(B*Delta)*B'*Ai,  -Ai*(B*Delta);
                       -(Delta*B')*Ai       ,    (Delta)     ]*Y;

    BetaOptimal1  =OptimalParametrs(0*rMp+0*rDp+1 : 1*rMp+0*rDp    );
    BetaOptimal3  =OptimalParametrs(1*rMp+1*rDp+1 : 2*rMp+1*rDp    );
    ThetaOptimal13=OptimalParametrs(2*rMp+2*rDp+1 : 2*rMp+2*rDp+rBp);

    EtaOptimal1   =0;
    EtaOptimal3   =0;

elseif(rM>0 && rD>-1 && rB>0) % mono+damper+bi

    i=1;
    A((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )    =inv([ (1+Gamma)*W1*(QQ_1'*QQ_1) , W1*(QQ_1'*DD_1);  W1*(DD_1'*QQ_1) ,(1+Gamma)*W1*(DD_1'*DD_1)]+Landa*eye(rMp+rDp));
    Amain((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )=([ (1+Gamma)*W1*(QQ_1'*QQ_1) , W1*(QQ_1'*DD_1);  W1*(DD_1'*QQ_1) ,(1+Gamma)*W1*(DD_1'*DD_1)]+Landa*eye(rMp+rDp));
    B((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , 1:rBp )=[W1*(QQ_1'*QH_13); W1*(DD_1'*QH_13)];
    Y((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , 1 )=W1*([QQ_1';DD_1']*Ur1(1:SampleRate:end));

    i=2;
    A((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )    =inv([ (1+Gamma)*W3*(QQ_3'*QQ_3) , W3*(QQ_3'*DD_3);  W3*(DD_3'*QQ_3) , (1+Gamma)*W3*(DD_3'*DD_3)]+Landa*eye(rMp+rDp));
    Amain((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , (i-1)*rMp+(i-1)*rDp+1:(i-0)*rMp+(i-0)*rDp )=([ (1+Gamma)*W3*(QQ_3'*QQ_3) , W3*(QQ_3'*DD_3);  W3*(DD_3'*QQ_3) , (1+Gamma)*W3*(DD_3'*DD_3)]+Landa*eye(rMp+rDp));
    B((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , 1:rBp )=[W3*(QQ_3'*QH_13); W3*(DD_3'*QH_13)];
    Y((i-1)*rMp+(i-1)*rDp+1 :(i-0)*rMp+(i-0)*rDp , 1 )=W3*([QQ_3';DD_3']*Ur3(1:SampleRate:end));

    D=(W1+W3)*(1+Gamma)*(QH_13'*QH_13)+Landa*eye(rBp);
    Y((2)*rMp+(2)*rDp+1 :(2)*rMp+(2)*rDp+rBp , 1 )=W1*(QH_13'*Ur1(1:SampleRate:end))+W3*(QH_13'*Ur3(1:SampleRate:end));
    
    Ai=A;
    Delta=inv(D-B'*Ai*B);
    OptimalParametrs=[ Ai+Ai*(B*Delta)*B'*Ai,  -Ai*(B*Delta);
                       -(Delta*B')*Ai       ,    (Delta)     ]*Y;

    BetaOptimal1  =OptimalParametrs(0*rMp+0*rDp+1 : 1*rMp+0*rDp    );
    BetaOptimal3  =OptimalParametrs(1*rMp+1*rDp+1 : 2*rMp+1*rDp    );
    ThetaOptimal13=OptimalParametrs(2*rMp+2*rDp+1 : 2*rMp+2*rDp+rBp);

    EtaOptimal1   =OptimalParametrs(1*rMp+0*rDp+1 : 1*rMp+1*rDp    );
    EtaOptimal1 =diag(1-sign(EtaOptimal1))/2*EtaOptimal1;
    EtaOptimal3   =OptimalParametrs(2*rMp+1*rDp+1 : 2*rMp+2*rDp    );
    EtaOptimal3 =diag(1-sign(EtaOptimal3))/2*EtaOptimal3;
    
    
    if(~estimate)
        H=[Amain,B;B',D];    
        options = optimoptions('quadprog',...
            'Algorithm','interior-point-convex','Display','off','MaxIter',100,...
            'TolFun',.00001);
        Xstar=quadprog(H,-Y,[],[],[],[],-inf*ones(2*rMp+2*rDp+rBp,1),[inf*ones(rMp,1);zeros(rDp,1);inf*ones(rMp,1);zeros(rDp,1);inf*ones(rBp,1)],OptimalParametrs,options);
    
        BetaOptimal1  =Xstar(0*rMp+0*rDp+1 : 1*rMp+0*rDp    );
        BetaOptimal3  =Xstar(1*rMp+1*rDp+1 : 2*rMp+1*rDp    );
        ThetaOptimal13=Xstar(2*rMp+2*rDp+1 : 2*rMp+2*rDp+rBp);

        EtaOptimal1   =Xstar(1*rMp+0*rDp+1 : 1*rMp+1*rDp    );
        EtaOptimal3   =Xstar(2*rMp+1*rDp+1 : 2*rMp+2*rDp    );
    end
    

end

% Torque
TorqueMonoSpring1 = QQc_1 *BetaOptimal1  ; 
TorqueMonoSpring3 = QQc_3 *BetaOptimal3  ; 
TorqueDamper1     = DDc_1 *EtaOptimal1   ; 
TorqueDamper3     = DDc_3 *EtaOptimal3   ; 
TorqueBiSpring13  = QHc_13*ThetaOptimal13; 


BetaOptimal =[BetaOptimal1;BetaOptimal3];
EtaOptimal  =[EtaOptimal1 ;EtaOptimal3 ];
ThetaOptimal= ThetaOptimal13;

% Cost
CostActuation=0;
CostD2Q=0;
CostParaReg=0;


Uactive1=(Ur1 -TorqueMonoSpring1 -TorqueDamper1 -TorqueBiSpring13 );
Uactive3=(Ur3 -TorqueMonoSpring3 -TorqueDamper3 -TorqueBiSpring13 );

CostActuation=1/2*sum(( W1*(Uactive1.*Uactive1) + W3*(Uactive3.*Uactive3)).*[diff(Time) ;0]);
% CostRequired=1/2*sum(( W1*(Ur1.*Ur1) + W3*(Ur3.*Ur3)).*[diff(Time) ;0]);
  
if(Show)
    warning off
    
    figure(20)
        clf
        ap=get(gca,'position');
        set(gcf,'position',[6 357 1266  195])
        % Required
        ax1=subplot(1,2,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(rad2deg(Q1),Ur1,'linewidth',2)
%         title('Hip joint')
        xlabel('\boldmath$q_1$ (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_1$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg(Q1),Uactive1,'linewidth',2,'color','r','linestyle','-.')
        hold off
        legend('u_r_1','u_a_1','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                  'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(Q1(1:150:SwitchIndex)),Ur1(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(Q1(SwitchIndex+1:75:end)),Ur1(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('Stance','Swing');
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1]);
        drawnow
        set(ax2,'visible','off','color','none',...
                 'box','off',...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        
        ax1=subplot(1,2,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(rad2deg( Q3),Ur3,'linewidth',2);
%         title('Knee joint')
        xlabel('\boldmath$q_3$ (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_3$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg(Q3),Uactive3,'linewidth',2,'color','r','linestyle','-.')
        hold off
        legend('u_r_3','u_a_3','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(Q3(1:150:SwitchIndex)),Ur3(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(Q3(SwitchIndex+1:75:end)),Ur3(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('Stance','Swing');
        drawnow
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1],...
                'fontweight','bold');
        set(ax2,'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[],...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        
    
    figure(21)
        % Mono-articular Spring
        % Mono-articular Spring
        % Mono-articular Spring
        ax1=subplot(3,2,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
        plot(rad2deg(Q1),TorqueMonoSpring1,'linewidth',3,'color','g')
        xlabel('\boldmath$q_1$ (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_1$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg( Q1),Ur1-TorqueBiSpring13-TorqueDamper1,'linewidth',2,'linestyle','-.')
        hold off
        legend('u_m_1','u_r_1-u_b_1-u_d_1','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(Q1(1:150:SwitchIndex)),Ur1(1:150:SwitchIndex)-TorqueBiSpring13(1:150:SwitchIndex)-TorqueDamper1(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(Q1(SwitchIndex+1:75:end)),Ur1(SwitchIndex+1:75:end)-TorqueBiSpring13(SwitchIndex+1:75:end)-TorqueDamper1(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('Stance','Swing');
        drawnow
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1]);
        set(ax2,'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[],...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        %
        ax1=subplot(3,2,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(rad2deg(Q3),TorqueMonoSpring3,'linewidth',3,'color','g')
        xlabel('\boldmath$q_3$ (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_3$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg(Q3),Ur3-TorqueDamper3-TorqueBiSpring13,'linewidth',2,'linestyle','-.')
        hold off
        legend('u_m_3','u_r_3-u_d_3-u_b_1_3','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(Q3(1:150:SwitchIndex)),Ur3(1:150:SwitchIndex)-TorqueBiSpring13(1:150:SwitchIndex)-TorqueDamper3(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(Q3(SwitchIndex+1:75:end)),Ur3(SwitchIndex+1:75:end)-TorqueBiSpring13(SwitchIndex+1:75:end)-TorqueDamper3(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('Stance','Swing');
        drawnow
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1]);
        set(ax2,'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[],...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        % Damper
        % Damper
        % Damper
        ax1=subplot(3,2,3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
        plot(rad2deg(DtQ1),TorqueDamper1,'linewidth',3,'Color',[1 0.7 .4])
        xlabel('\boldmath$\dot{q}_1$ (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_1$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg(DtQ1),Ur1-TorqueBiSpring13-TorqueMonoSpring1,'linewidth',2,'linestyle','-.')
        hold off
        legend('u_d_1','u_r_1-u_b_1-u_m_1','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(DtQ1(1:150:SwitchIndex)),Ur1(1:150:SwitchIndex)-TorqueBiSpring13(1:150:SwitchIndex)-TorqueMonoSpring1(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(DtQ1(SwitchIndex+1:75:end)),Ur1(SwitchIndex+1:75:end)-TorqueBiSpring13(SwitchIndex+1:75:end)-TorqueMonoSpring1(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('Stance','Swing');
        drawnow
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1],...
                'fontweight','bold');
        set(ax2,'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[],...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        %
        ax1=subplot(3,2,4,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(rad2deg(DtQ3),TorqueDamper3,'linewidth',3,'Color',[1 0.7 .4])
        xlabel('\boldmath$\dot{q}_3$ (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_3$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg(DtQ3),Ur3-TorqueMonoSpring3-TorqueBiSpring13,'linewidth',2,'linestyle','-.')
        hold off
        legend('u_d_3','u_r_3-u_m_3-u_b_1_3','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(DtQ3(1:150:SwitchIndex)),Ur3(1:150:SwitchIndex)-TorqueBiSpring13(1:150:SwitchIndex)-TorqueMonoSpring3(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(DtQ3(SwitchIndex+1:75:end)),Ur3(SwitchIndex+1:75:end)-TorqueBiSpring13(SwitchIndex+1:75:end)-TorqueMonoSpring3(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('stance','Swing');
        drawnow
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1],...
                'fontweight','bold');
        set(ax2,'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[],...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        % Bi-articular Spring
        % Bi-articular Spring
        % Bi-articular Spring
        ax1=subplot(325);   
        sp1=get(ax1,'position');
        set(ax1,'position',[sp1(1)+.5*(ap(3)-sp1(3)),sp1(2:end)]); 
        plot(rad2deg(Q1+Q3),2*TorqueBiSpring13,'linewidth',3,'Color',[0.75 0 0.75])
        xlabel('\boldmath$q_1+q_3$ (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        ylabel('\boldmath$u_1+u_3$ (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10','interpret','latex');
        grid on
        set(gca,'YMinorGrid','on')
        hold on
        plot(rad2deg(Q1+Q3),Ur1+Ur3-TorqueMonoSpring1-TorqueMonoSpring3-TorqueDamper1-TorqueDamper3,'linewidth',2,'linestyle','-.')
        hold off
        legend('2\times u_b_1_3','u_r_1+u_r_3-u_m_1-u_m_3-u_d_1-u_d_3','Orientation','horizontal')
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ax2=axes('position',get(ax1,'position'),'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[]);
        plot(rad2deg(Q1(1:150:SwitchIndex)+Q3(1:150:SwitchIndex)),Ur1(1:150:SwitchIndex)+Ur3(1:150:SwitchIndex)-TorqueMonoSpring1(1:150:SwitchIndex)-TorqueMonoSpring3(1:150:SwitchIndex)-TorqueDamper1(1:150:SwitchIndex)-TorqueDamper3(1:150:SwitchIndex),'linestyle','none','marker','o','markerfacecolor','b')
        hold on
        plot(rad2deg(Q1(SwitchIndex+1:75:end)+Q3(SwitchIndex+1:75:end)),Ur1(SwitchIndex+1:75:end)+Ur3(SwitchIndex+1:75:end)-TorqueMonoSpring1(SwitchIndex+1:75:end)-TorqueMonoSpring3(SwitchIndex+1:75:end)-TorqueDamper1(SwitchIndex+1:75:end)-TorqueDamper3(SwitchIndex+1:75:end),'linestyle','none','marker','^','markerfacecolor','b')        
        L2=legend('Stance','Swing');
        drawnow
        set(L2,'Orientation','horizontal','Location','SouthWest','Color',[1 1 1],...
                'fontweight','bold');
        set(ax2,'visible','off','color','none',...
                 'box','off','xtick',[],'ytick',[],...
                 'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),...
                  'position',get(ax1,'position'));        
        


    figure(22)
        clf
        subplot(2,2,1)
        plot(rad2deg(DtQ1),TorqueDamper1,'linewidth',2)
        hold all
        plot(rad2deg(DtQ1),Ur1-TorqueMonoSpring1-TorqueBiSpring13,'-.','linewidth',2)
        hold off
        grid on
        xlabel('${\dot q_1}$ (deg/s)', 'interpreter','latex','fontsize',14,'FontName','mwa_cmb10');
        ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')    
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        legend('u_d_1','u_r_1-u_m_1-u_b_1_3','Orientation','horizontal')
        subplot(2,2,2)
        plot(rad2deg(DtQ1),DtQ1.*TorqueDamper1,'linewidth',2)
        grid on
        xlabel('${\dot q_1}$ (deg/s)', 'interpreter','latex','fontsize',14,'FontName','mwa_cmb10');
        ylabel('${u_1 * \dot q_1}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');                
        set(gca,'YMinorGrid','on')    
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%         legend('u_d_1','u_r,'Orientation','horizontal')
        title('power')
        subplot(2,2,3)
        plot(rad2deg(DtQ3),TorqueDamper3,'linewidth',2)
        hold all
        plot(rad2deg(DtQ3),Ur3-TorqueMonoSpring3-TorqueBiSpring13,'-.','linewidth',2)
        hold off
        xlabel('${\dot q_3}$ (deg/s)', 'interpreter','latex','fontsize',14,'FontName','mwa_cmb10');
        ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')    
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        legend('u_d_3','u_r_3-u_m_3-u_b_1_3','Orientation','horizontal')
        grid on
        subplot(2,2,4)
        plot(rad2deg(DtQ3),DtQ3.*TorqueDamper3,'linewidth',2)
        grid on
        xlabel('${\dot q_3}$ (deg/s)', 'interpreter','latex','fontsize',14,'FontName','mwa_cmb10');
        ylabel('${u_3 * \dot q_3}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');                
        set(gca,'YMinorGrid','on')    
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%         legend('u_d_1','u_r,'Orientation','horizontal')

    warning on
end

end