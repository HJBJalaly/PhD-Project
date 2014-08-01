function DQ=SirDynSample(t,Q,g,L,m,KKmat)
    
q1=Q(1);
q2=Q(2);
q3=Q(3);

D1q1=Q(4);
D1q2=Q(5);
D1q3=Q(6);
D1q=[D1q1 D1q2 D1q3]';

mL1=m;
mL2=m;
mL3=m;

LL1=L;
LL2=L;
LL3=L;


% Dynamic
MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
    mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
    mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];


CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
      mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
      mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];


GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
     g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
     g*mL3*LL3*cos(q1+q2+q3)/2];

Gamma=-KKmat.*Q(1:3);


% Apply Torqe

D2q= MM^-1*(-CC*D1q-GG + Gamma);
De=(Gamma'*D1q);
DQ=[D1q;D2q;De];

end