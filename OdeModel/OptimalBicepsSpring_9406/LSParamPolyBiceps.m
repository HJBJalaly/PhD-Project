function  [BetaOptimal,ThetaOptimal,Cost,TorquePassiveOptimal,TorqueBicepsOptimal]=...
                    LSParamPolyBiceps(Qq,Qhatq,TorqueDesire,nn,rU,rB,Landa,SampleRate)

                
                
if (rB>0)                
    QQ={};
    QQhat={};
    D2Q={};
    D2Qhat={};
    

    % for 1 to n-1
    for kk=1:nn-1    
        Q_i=[];
        Qhat_i=[];
        D2Q_i=[];
        D2Qhat_i=[];
        for i=1:SampleRate:size(Qq,2)
            Q_i(i,:) = Qq(kk,i).^(rU:-1:0)';
            Qhat_i(i,:) = Qhatq(kk,i).^(rB:-1:0)';
            
            D2Q_i(i,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
            D2Qhat_i(i,:) = ([Qhatq(kk,i).^(rB-2:-1:0) 0 0].*(rB:-1:0).*(rB-1:-1:-1))';
            
        end
        QQ{kk}=Q_i;
        QQhat{kk}=Qhat_i;
        D2Q{kk}=D2Q_i;
        D2Qhat{kk}=D2Qhat_i;
    end
    % for n
    Q_i=[];
    D2Q_i=[];
    for i=1:SampleRate:size(Qq,2)
        Q_i(end+1,:) = Qq(nn,i).^(rU:-1:0)';
        D2Q_i(i,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
    end
    QQ{nn}=Q_i;
    D2Q{nn}=D2Q_i;


    A=[];
    B=[];
    D=[];
    Y=[];

    rUp=rU+1;
    rBp=rB+1;
    % for i=1 to n-2
    for i=1:nn-2
        A((i-1)*rUp+1:(i)*rUp,(i-1)*rUp+1:(i)*rUp)=inv(QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*ones(size(QQ{i}'*QQ{i})));

        B((i-1)*rUp+1:(i)*rUp,(i-1)*rBp+1:(i)*rBp)=QQ{i}'*QQhat{i};
        B((i)*rUp+1:(i+1)*rUp,(i-1)*rBp+1:(i)*rBp)=QQ{i+1}'*QQhat{i};

        D((i-1)*rBp+1:(i)*rBp,(i-1)*rBp+1:(i)*rBp)=2*QQhat{i}'*QQhat{i}+Landa(4)*D2Qhat{i}'*D2Qhat{i}+Landa(3)*eye(size(QQhat{i}'*QQhat{i}));
        D((i-1)*rBp+1:(i)*rBp,(i)*rBp+1:(i+1)*rBp)=QQhat{i}'*QQhat{i+1};
        D((i)*rBp+1:(i+1)*rBp,(i-1)*rBp+1:(i)*rBp)=QQhat{i+1}'*QQhat{i};

        Y((i-1)*rUp+1:(i)*rUp,1)=QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
        Y(nn*rUp+(i-1)*rBp+1:(nn*rUp)+(i)*rBp,1)=QQhat{i}'*(TorqueDesire(i,1:SampleRate:end)'+TorqueDesire(i+1,1:SampleRate:end)');
    end
    % for i=n-1
        i=nn-1;
        A((i-1)*rUp+1:(i)*rUp,(i-1)*rUp+1:(i)*rUp)=inv(QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*eye(size(QQ{i}'*QQ{i})));

        B((i-1)*rUp+1:(i)*rUp,(i-1)*rBp+1:(i)*rBp)=QQ{i}'*QQhat{i};
        B((i)*rUp+1:(i+1)*rUp,(i-1)*rBp+1:(i)*rBp)=QQ{i+1}'*QQhat{i};

        D((i-1)*rBp+1:(i)*rBp,(i-1)*rBp+1:(i)*rBp)=2*QQhat{i}'*QQhat{i}+Landa(4)*D2Qhat{i}'*D2Qhat{i}+Landa(3)*eye(size(QQhat{i}'*QQhat{i}));

        Y((i-1)*rUp+1:(i)*rUp,1)=QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
        Y(nn*rUp+(i-1)*rBp+1:(nn*rUp)+(i)*rBp,1)=QQhat{i}'*(TorqueDesire(i,1:SampleRate:end)'+TorqueDesire(i+1,1:SampleRate:end)');
    % for i=n
        i=nn;
        A((i-1)*rUp+1:(i)*rUp,(i-1)*rUp+1:(i)*rUp)=inv(QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*eye(size(QQ{i}'*QQ{i})));

        Y((i-1)*rUp+1:(i)*rUp,1)=QQ{i}'*TorqueDesire(i,1:SampleRate:end)';


    
    Ai=(A);
    Delta=inv(D-B'*Ai*B);

    OptimalSpring=[Ai+Ai*(B*Delta)*B'*Ai,  -Ai*(B*Delta);
                 -(Delta*B')*Ai       ,   (Delta)        ]*Y;

    BetaOptimal=OptimalSpring(1:nn*(rU+1));
    ThetaOptimal=OptimalSpring(nn*(rU+1)+1:end);

    % Cost
    Cost=0;

    i=1;
    Cost=Cost + ...
          ( 1/2*(TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rU+1))-QQhat{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)) )'*...
                (TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhat{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)) )+...
                Landa(2)* BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Q{i}'*D2Q{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))+...
                Landa(4)* ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))'*(D2Qhat{i}'*D2Qhat{i})*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)));
    for i=2:nn-1

        Cost=Cost + ...
          ( 1/2*(TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhat{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))-QQhat{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp)))'*...
                (TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhat{i}*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))-QQhat{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp)))+...
                Landa(2)* BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Q{i}'*D2Q{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))+...
                Landa(4)* ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp))'*(D2Qhat{i}'*D2Qhat{i})*ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)));
    end
    i=nn;
    Cost=Cost + ...
          ( 1/2*(TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhat{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp)))'*...
                (TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))-QQhat{i-1}*ThetaOptimal((i-2)*(rBp)+1:(i-1)*(rBp))))+...
                Landa(2)* BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Q{i}'*D2Q{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp));
    % Torque

    for i=1:nn-1
        TorquePassiveOptimal(i,:)=polyval(BetaOptimal((i-1)*(rUp)+1:i*(rUp)),Qq(i,:));

        TorqueBicepsOptimal(i,:)=polyval(ThetaOptimal((i-1)*(rBp)+1:(i)*(rBp)),Qhatq(i,:));
    end
    i=nn;
    TorquePassiveOptimal(i,:)=polyval(BetaOptimal((i-1)*(rUp)+1:i*(rUp)),Qq(i,:));
    
else
%%    
    QQ={};
    D2Q={};
    
    % w
    % for 1 to n-1
    for kk=1:nn
        Q_i=[];
        D2Q_i=[];

        for i=1:SampleRate:size(Qq,2)
            Q_i(end+1,:) = Qq(kk,i).^(rU:-1:0)';
            D2Q_i(i,:) = ([Qq(kk,i).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
        end
        QQ{kk}=Q_i;
        D2Q{kk}=D2Q_i;
    end
    

    % 
    A=[];
    
    rUp=rU+1;
    % for i=1 to n-2
    for i=1:nn
        A((i-1)*rUp+1:(i)*rUp,(i-1)*rUp+1:(i)*rUp)=(QQ{i}'*QQ{i}+Landa(2)*D2Q{i}'*D2Q{i}+Landa(1)*ones(size(QQ{i}'*QQ{i})));


        Y((i-1)*rUp+1:(i)*rUp,1)=QQ{i}'*TorqueDesire(i,1:SampleRate:end)';
    end
    
    OptimalSpring=A\Y(1:nn*(rUp));
    BetaOptimal=OptimalSpring(1:nn*(rUp));
    ThetaOptimal=[];

    % Cost



    Cost=0;

    for i=1:nn

        Cost=Cost + ...
          ( 1/2*(TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp)))'*...
                (TorqueDesire(i,:)'- QQ{i}*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))))+...
                Landa(2)* BetaOptimal((i-1)*(rUp)+1:(i)*(rUp))'*(D2Q{i}'*D2Q{i})*BetaOptimal((i-1)*(rUp)+1:(i)*(rUp));
    end
    
    % Torque

    for i=1:nn
        TorquePassiveOptimal(i,:)=polyval(BetaOptimal((i-1)*(rUp)+1:i*(rUp)),Qq(i,:));
        
    end
    TorqueBicepsOptimal=zeros(nn-1,size(Qq,2));
end


%%
% 
% if(rU>0)
%     if(rB>0)
%         Temp=(eye(length(PQ))- PQ*((PQ'*PQ)\PQ'));
%     else
%         Temp=eye(length(PQ));
%     end
%     Beta=((QQ'*Temp*QQ))\( QQ'*Temp*Torque(1:SampleRate:end));
% else
%     Beta=[0 0];
% end
% 
% 
% if(rB>0)
%     if(rU>0)
% %         Gamma=(eye(length(PQ))- PQ*((PQ'*PQ)\PQ'));
%         Theta=(PQ'*PQ)\( PQ'*Torque(1:SampleRate:end)-PQ'*QQ*Beta );
%     else
%         Theta=((PQ'*PQ))\( PQ'*Torque(1:SampleRate:end));
%     end
% else
%     Theta=[0 0];
% end


end