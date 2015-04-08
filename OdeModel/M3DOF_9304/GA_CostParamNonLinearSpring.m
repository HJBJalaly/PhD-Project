function Cost = GA_CostParamNonLinearSpring(Param,ThetaStep,ThetaS,tau)
% Param
% ThetaStep : in rad
% Theta: in rad with constant diff
% tau: in N.m

% if(any(isinf( Param)))
%     1;
% end

k=Param(1);
R=Param(2);
q0=Param(3);

sum=0;
for i=1:length(ThetaS)
    J(i)=tau(i)/sqrt(2*k*sum+(k*q0)^2);
    sum=sum+tau(i)*(ThetaStep);
end

DJ=differential(J,ThetaS,(ThetaStep));
r= sqrt( J.^2 + (DJ.^2 .* (R^2 - J.^2))./((DJ+sqrt(R^2-J.^2)).^2) );
ThetaR=-ThetaS+acos((J.^2+sign(DJ).*sqrt((r.^2-J.^2).*(R^2-J.^2))  )./(R*r));

PenTR=max(abs(diff(sign(diff(ThetaR)))));
Cost= k/100 + 5*q0 + 2*R + 100*(max(abs(diff(r)))) + 4*(R-mean(r)) + 1000*PenTR + 1500*(R<1.5*max(r))  + ...
            1e6*(~(isreal(r)))+ 5e6*(~(isreal(R)))+ 5e5*(~(isreal(k)))+ 2e5*(~(isreal(q0)))+ 8e5*any(Param<0);

% if(isinf( Cost))
%     1;
% end

end

