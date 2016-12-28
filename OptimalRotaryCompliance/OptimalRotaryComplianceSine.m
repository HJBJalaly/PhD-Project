function  [BetaOptimal,CostActuation,CostD2Q,Up]=...
                    OptimalRotaryComplianceSine(Qq,Ur,rp,gamma,omega,delta,SampleRate,scale)

    Qq=Qq/scale;
    PHI=[];
    DD=zeros(rp+1);
    PHIc=[];
    DDc=zeros(rp+1);

    for i=1:SampleRate:length(Qq)
        PHI(:,end+1) = Qq(i).^(rp:-1:0)';
        D = ([Qq(i).^(rp-2:-1:0) 0 0].*(rp:-1:0).*(rp-1:-1:-1))';
        DD=DD+D*D';
    end
    for i=1:length(Qq)
        PHIc(:,end+1) = Qq(i).^(rp:-1:0)';
        D = ([Qq(i).^(rp-2:-1:0) 0 0].*(rp:-1:0).*(rp-1:-1:-1))';
        DDc=DDc+D*D';
    end
        
    Sigma=PHI*PHI'+omega*gamma*DD+delta*omega*eye(rp+1);
    
    A=[];
    for i=rp:-1:0
        A(end+1,1)=1/(i+1)*(2*pi/scale)^(i+1);
    end
    
    Ur_s=Ur(1:SampleRate:end);
    BetaOptimal=Sigma\(...  
                       PHI*Ur_s- ((A'*Sigma^(-1)*A)^(-1))*(A'*Sigma^(-1)*PHI*Ur_s*A) ...
                        );
    
    % Cost

    CostActuation= 1/2*(Ur- PHIc'*BetaOptimal)'*(Ur- PHIc'*BetaOptimal);
    CostD2Q= BetaOptimal'*DD*BetaOptimal;
    
    % Torque

    Up=polyval(BetaOptimal,Qq);
        
end
