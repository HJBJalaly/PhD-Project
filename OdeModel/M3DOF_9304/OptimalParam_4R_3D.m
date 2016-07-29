function  [BetaOptimal,ThetaOptimal,CostActuation,CostD2Q,CostParaReg,TorqueMonoOptimal,TorqueBicepsOptimal]=...
                    OptimalParam_4R_3D(Qq,Qhatq,TorqueDesire,nn,rU,rB,Landa,Weight,SampleRate)

                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rU>0 : Filling Trajectory Matrices
if (rB>0)                
    QQ={};
    QQhat={};
    D2Q={};
    D2Qhat={};
    QQc={};
    QQhatc={};
    D2Qc={};
    D2Qhatc={};
    

    % for 1 to n-1
    for kk=1:nn-1    
        Q_i=[];
        Qhat_i=[];
        D2Q_i=[];
        D2Qhat_i=[];
        Qc_i=[];
        Qhatc_i=[];
        D2Qc_i=[];
        D2Qhatc_i=[];
        for i=1:SampleRate:size(Qq,2)
            Q_i(end+1,:) = Qq(kk,i).^(rU:-1:0)';
            Qhat_i(end+1,:) = Qhatq(kk,i).^(rB:-1:0)';
            
            D2Q_i(end+1,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
            D2Qhat_i(end+1,:) = ([Qhatq(kk,i).^(rB-2:-1:0) 0 0].*(rB:-1:0).*(rB-1:-1:-1))';
            
        end
        for i=1:size(Qq,2)
            Qc_i(end+1,:) = Qq(kk,i).^(rU:-1:0)';
            Qhatc_i(end+1,:) = Qhatq(kk,i).^(rB:-1:0)';
            
            D2Qc_i(end+1,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
            D2Qhatc_i(end+1,:) = ([Qhatq(kk,i).^(rB-2:-1:0) 0 0].*(rB:-1:0).*(rB-1:-1:-1))';
            
        end
        QQ{kk}=Q_i;
        QQhat{kk}=Qhat_i;
        D2Q{kk}=D2Q_i;
        D2Qhat{kk}=D2Qhat_i;
        
        QQc{kk}=Qc_i;
        QQhatc{kk}=Qhatc_i;
        D2Qc{kk}=D2Qc_i;
        D2Qhatc{kk}=D2Qhatc_i;
    end
    % for n
    Q_i=[];
    D2Q_i=[];
    Qc_i=[];
    D2Qc_i=[];
    for i=1:SampleRate:size(Qq,2)
        Q_i(end+1,:) = Qq(nn,i).^(rU:-1:0)';
        D2Q_i(end+1,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
    end
    for i=1:size(Qq,2)
        Qc_i(end+1,:) = Qq(nn,i).^(rU:-1:0)';
        D2Qc_i(end+1,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
    end
    QQ{nn}=Q_i;
    D2Q{nn}=D2Q_i;
    QQc{nn}=Qc_i;
    D2Qc{nn}=D2Qc_i;


    
    rUp=rU+1;
    rBp=rB+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% rB>0 : Optimize the compliance of first joint
    
    A=[];
    Y=[];

    i=1;
    A((i-1)*rUp+1:(i)*rUp,(i-1)*rUp+1:(i)*rUp)=(Weight(i)*QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*((QQ{i}'*QQ{i})));
    Y((i-1)*rUp+1:(i)*rUp,1)=Weight(i)*QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
    
    OptimalSpring_FirstJoint=A\Y;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% rB>0 : Optimize the compliance of first joint
    A=[];
    B=[];
    D=[];
    Y=[];

    NumberOfLastJoints=nn-1;
    
    % for i=1 to n-2
    for i=2:nn-2
        im=i-1;
        A((im-1)*rUp+1:(im)*rUp,(im-1)*rUp+1:(im)*rUp)=inv(Weight(i)*QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*((QQ{i}'*QQ{i})));

        B((im-1)*rUp+1:(im)*rUp,(im-1)*rBp+1:(im)*rBp)=Weight(i)*QQ{i}'*QQhat{i};
        B((im)*rUp+1:(im+1)*rUp,(im-1)*rBp+1:(im)*rBp)=Weight(i+1)*QQ{i+1}'*QQhat{i};

        D((im-1)*rBp+1:(im)*rBp,(im-1)*rBp+1:(im)*rBp)=(Weight(i)+Weight(i+1))*QQhat{i}'*QQhat{i}+Landa(4)*D2Qhat{i}'*D2Qhat{i}+Landa(3)*((QQhat{i}'*QQhat{i}));
        D((im-1)*rBp+1:(im)*rBp,(im)*rBp+1:(im+1)*rBp)=Weight(i+1)*QQhat{i}'*QQhat{i+1};
        D((im)*rBp+1:(im+1)*rBp,(im-1)*rBp+1:(im)*rBp)=Weight(i+1)*QQhat{i+1}'*QQhat{i};

        Y((im-1)*rUp+1:(im)*rUp,1)=Weight(i)*QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
        Y((NumberOfLastJoints)*rUp+(im-1)*rBp+1:((NumberOfLastJoints)*rUp)+(im)*rBp,1)=QQhat{i}'*(Weight(i)*TorqueDesire(i,1:SampleRate:end)'+Weight(i+1)*TorqueDesire(i+1,1:SampleRate:end)');
    end
    % for i=n-1
    i=nn-1;
    im=i-1;
        A((im-1)*rUp+1:(im)*rUp,(im-1)*rUp+1:(im)*rUp)=inv(Weight(i)*QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*((QQ{i}'*QQ{i})));

        B((im-1)*rUp+1:(im)*rUp,(im-1)*rBp+1:(im)*rBp)=Weight(i)*QQ{i}'*QQhat{i};
        B((im)*rUp+1:(im+1)*rUp,(im-1)*rBp+1:(im)*rBp)=Weight(i+1)*QQ{i+1}'*QQhat{i};

        D((im-1)*rBp+1:(im)*rBp,(im-1)*rBp+1:(im)*rBp)=(Weight(i)+Weight(i+1))*QQhat{i}'*QQhat{i}+Landa(4)*D2Qhat{i}'*D2Qhat{i}+Landa(3)*((QQhat{i}'*QQhat{i}));

        Y((im-1)*rUp+1:(im)*rUp,1)=Weight(i)*QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
        Y((NumberOfLastJoints)*rUp+(im-1)*rBp+1:((NumberOfLastJoints)*rUp)+(im)*rBp,1)=QQhat{i}'*(Weight(i)*TorqueDesire(i,1:SampleRate:end)'+Weight(i+1)*TorqueDesire(i+1,1:SampleRate:end)');
    % for i=n
    i=nn;
    im=i-1;
        A((im-1)*rUp+1:(im)*rUp,(im-1)*rUp+1:(im)*rUp)=inv(Weight(i)*QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*((QQ{i}'*QQ{i})));

        Y((im-1)*rUp+1:(im)*rUp,1)=Weight(i)*QQ{i}'*TorqueDesire(i,1:SampleRate:end)';


    
    Ai=(A);
    Delta=inv(D-B'*Ai*B);

    OptimalSpring_Last3Joint=[Ai+Ai*(B*Delta)*B'*Ai,  -Ai*(B*Delta);
                              -(Delta*B')*Ai       ,   (Delta)      ]*Y;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    BetaOptimal=[OptimalSpring_FirstJoint ;OptimalSpring_Last3Joint(1:(NumberOfLastJoints)*(rUp))];
    ThetaOptimal=[zeros(rBp,1); OptimalSpring_Last3Joint((NumberOfLastJoints)*(rUp)+1:end)];

    
    
    
    % Cost
    CostActuation=0;
    CostD2Q=0;
    CostParaReg=0;


    i=1;
    CostActuation=CostActuation + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rU+1))-QQhatc{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)) )'*...
                    (TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhatc{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)) ));
                    
    CostD2Q=CostD2Q + ... 
          ([ BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Qc{i}'*D2Qc{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp)),
                       ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))'*(D2Qhatc{i}'*D2Qhatc{i})*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))]);
    
    for i=2:nn-1

        CostActuation=CostActuation + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhatc{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))-QQhatc{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp)))'*...
                    (TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhatc{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))-QQhatc{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp))));
        CostD2Q=CostD2Q + ...
          Weight(i)*([ BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Qc{i}'*D2Qc{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp)),
                       ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))'*(D2Qhatc{i}'*D2Qhatc{i})*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))]);
    end
    i=nn;
    CostActuation=CostActuation + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhatc{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp)))'*...
                    (TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhatc{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp))));
    CostD2Q=CostD2Q + ...
          Weight(i)*([ BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Qc{i}'*D2Qc{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp)),
                        0]);
    
    CostParaReg=[ BetaOptimal'*BetaOptimal,
                  ThetaOptimal'*ThetaOptimal ];
                
    % Torque
    
    for i=1:nn-1
        TorqueMonoOptimal(i,:)=polyval(BetaOptimal((i-1)*(rUp)+1:i*(rUp)),Qq(i,:));

        TorqueBicepsOptimal(i,:)=polyval(ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)),Qhatq(i,:));
    end
    i=nn;
    TorqueMonoOptimal(i,:)=polyval(BetaOptimal((i-1)*(rUp)+1:i*(rUp)),Qq(i,:));
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  %% rB=0 : Optimize the compliance of first joint
    
        
    
    QQ={};
    D2Q={};
    QQc={};
    QQhatc={};
    
    
    % w
    % for 1 to n-1
    for kk=1:nn
        Q_i=[];
        D2Q_i=[];
        Qc_i=[];
        D2Qc_i=[];
        

        for i=1:SampleRate:size(Qq,2)
            Q_i(end+1,:) = Qq(kk,i).^(rU:-1:0)';
            D2Q_i(end+1,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
        end
        for i=1:size(Qq,2)
            Qc_i(end+1,:) = Qq(kk,i).^(rU:-1:0)';
            D2Qc_i(end+1,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
        end
        QQ{kk}=Q_i;
        D2Q{kk}=D2Q_i;
        QQc{kk}=Qc_i;
        D2Qc{kk}=D2Qc_i;
    end
    

    % 
    A=[];
    
    rUp=rU+1;
    EYEU=[zeros(rU,rUp); zeros(1,rU),1 ];
    
    % for i=1 to n-2
    for i=1:nn
        A((i-1)*rUp+1:(i)*rUp,(i-1)*rUp+1:(i)*rUp)=(Weight(i)*QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*eye(size(QQ{i}'*QQ{i})));


        Y((i-1)*rUp+1:(i)*rUp,1)=Weight(i)*QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
    end
    
    OptimalSpring=A\Y(1:nn*(rUp));
    BetaOptimal=OptimalSpring(1:nn*(rUp));
    ThetaOptimal=[];

    % Cost



    CostActuation=0;
    CostD2Q=0;

    for i=1:nn

        CostActuation=CostActuation + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp)))'*...
                    (TorqueDesire(i,:)'- QQc{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))));
    
        CostD2Q=CostD2Q + ...
          Weight(i)*( [ Landa(2)* BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Qc{i}'*D2Qc{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp));
                        0]);
    end
    
    CostParaReg = [BetaOptimal'*BetaOptimal;  0 ];
    
    % Torque

    for i=1:nn
        TorqueMonoOptimal(i,:)=polyval(BetaOptimal((i-1)*(rUp)+1:i*(rUp)),Qq(i,:));
        
    end
    TorqueBicepsOptimal=zeros(nn-1,size(Qq,2));
end

end