function Cost=CF2_TorqueCost(Alpha,Time,Degree,Tres,Weight,Landa,QQ,B,g,mL1,mL2,mL3,LL1,LL2,LL3,MinSigValue)

% Torque=zeros(3,length(Time));
nn=Degree(1);
rQ=Degree(2);
rU=Degree(3);

Alpha_Q1=Alpha(1:(rQ+1));
Alpha_Q2=Alpha((rQ+1)+1:2*(rQ+1));
Alpha_Q3=Alpha(2*(rQ+1)+1:3*(rQ+1));


Alpha_D1Q1=Alpha_Q1(1:end-1).*(rQ:-1:1);
Alpha_D1Q2=Alpha_Q2(1:end-1).*(rQ:-1:1);
Alpha_D1Q3=Alpha_Q3(1:end-1).*(rQ:-1:1);

Alpha_D2Q1=Alpha_D1Q1(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q2=Alpha_D1Q2(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q3=Alpha_D1Q3(1:end-1).*(rQ-1:-1:1);

Q1=polyval(Alpha_Q1,Time);
Q2=polyval(Alpha_Q2,Time);
Q3=polyval(Alpha_Q3,Time);
QVal=[Q1;Q2;Q3];


D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);

D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);

for i=1:length(Time)
    q1=Q1(i);
    q2=Q2(i);
    q3=Q3(i);
    D1q1=D1Q1(i);
    D1q2=D1Q2(i);
    D1q3=D1Q3(i);
    D2q1=D2Q1(i);
    D2q2=D2Q2(i);
    D2q3=D2Q3(i);

    
    MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
        mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
        mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];

    CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
          mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
          mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];

    GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
         g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
         g*mL3*LL3*cos(q1+q2+q3)/2];


     TorqueDesire(:,i) = MM*[D2q1;D2q2;D2q3] + CC*[D1q1;D1q2;D1q3] + GG;
    
end



% Omega matrix
Omega = B * (B'*B)^-1 *diag(Weight) * (B'*B)^-1 * B';

% Integral Matrix
Iu=0;
Iq=zeros(nn*(rU+1),nn*(rU+1));
Iuq=zeros(1,(rU+1)*nn);
for tt=1:length(Time)
   
    QQ_conc=zeros((rU+1)*nn,nn);% \underline{\underline{\mathcal{Q}}}^{r_u}
    for joint=1:nn
        QQ_rU_Joint=QVal(joint,tt).^(rU:-1:0)';
        QQ_conc(1+(joint-1)*(rU+1):(joint)*(rU+1),joint)= QQ_rU_Joint;
    end

    Iu = Iu + TorqueDesire(:,tt)'*Omega*TorqueDesire(:,tt);

    Iq = Iq + QQ_conc*Omega*QQ_conc';

    Iuq= Iuq+ TorqueDesire(:,tt)'*Omega*QQ_conc';
   
end

Iu=Iu*Tres;
Iq=Iq*Tres;
Iuq=Iuq*Tres;


Idq_conc=zeros((rU+1)*nn,(rU+1)*nn);

for Joint=1:nn
    c_hat =(max(QVal(Joint,:)) - min(QVal(Joint,:)))* 2 / 3/pi;
    d_hat = min(QVal(Joint,:));

    for kk=1:rU+1
        Psi(kk,:) = [zeros(1,kk-1), (c_hat^(rU+1-(kk)))* poly(-d_hat/c_hat*ones(1,rU+1-(kk)))];
    end

    
    Idq = Weight(Joint)*c_hat* Psi*QQ*Psi';
    
    Idq_conc((Joint-1)*(rU+1)+1:(Joint)*(rU+1),(Joint-1)*(rU+1)+1:(Joint)*(rU+1)) = Idq;
end

% BetaOptimal=  Landa*(Landa*Iq + (1-Landa)*Idq_conc )^-1 *Iuq';
    
% IntU2=1/2*Iu+1/2*BetaOptimal'*Iq*BetaOptimal-Iuq*BetaOptimal;
% CostSlope=1/2*BetaOptimal'*Idq_conc*BetaOptimal;
% Cost = IntU2*Landa+(1-Landa)*CostSlope;
% Cost = 1/2* Landa*Iu - 1/2*Landa^2*Iuq*(Landa*Iq + (1-Landa)*Idq_conc )^-1 *Iuq';
SVDsol=SVDBlockInvertor((Landa*Iq+(1-Landa)*Idq_conc),nn,rU+1,MinSigValue);
Cost = 1/2* Landa*Iu - 1/2*Landa^2*Iuq*SVDsol*Iuq';
1;
end