function DQ=SirDyn2DoFPassive(t,Q,g,L,m,K1,K2)


q1=Q(1);
q2=Q(2);


D1q1=Q(3);
D1q2=Q(4);

D1q=[D1q1 D1q2 ]';

mL1=m;
mL2=m;

LL1=L;
LL2=L;


% Dynamic
MM=[mL2*LL1*LL2*cos(q2)+mL2*LL1^2+mL2*LL2^2/3+mL1*LL1^2/3, mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6;
    mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6, mL2*LL2^2/3];


CC= [-mL2*LL2*sin(q2)*D1q2*LL1/2, -mL2*LL2*sin(q2)*(D1q1+D1q2)*LL1/2; 
      mL2*LL2*LL1*D1q1*sin(q2)/2, 0];


GG= [g*(mL2*LL2*cos(q1+q2)+LL1*cos(q1)*mL1+2*LL1*cos(q1)*mL2)/2 + K1*q1;
     g*mL2*LL2*cos(q1+q2)/2 + K2*q2];

%  Control Law
Gamma=[0;0];
% Apply Torqe
D2q= MM^-1*(-CC*D1q-GG+Gamma);
De=abs(Gamma'*D1q);
DQ=[D1q;D2q;De];

end