function Dx=OdeSolverNonLinTRotExact(t,x,tau,Dtheta,D2theta,Time,K,l0,m)
                                     

Dx=zeros(2,1);

Tau     =interp1(Time,tau,t);
DTheta  =interp1(Time,Dtheta,t);
D2Theta =interp1(Time,D2theta,t);


Dx(1)=x(2);
Dx(2)=(2/3/m*(-K*(x(1)-l0)*x(2) - Tau* DTheta)       -  x(1)^2*DTheta*D2Theta - x(1)*x(2)*DTheta^2)    /x(2);

% if(Tau>=0)
%     Dx(2)=(2/3/m*(-K*(x(1)-l0)*x(2) + Tau* DTheta)       -  x(1)^2*DTheta*D2Theta - x(1)*x(2)*DTheta^2)    /x(2);
%  elseif(Tau==0)
%      Dx(1)=0;
%      Dx(2)=( 2/3/m*(-K*(x(1)-l0))       -  x(1)^2*DTheta*D2Theta/x(2) - x(1)*DTheta^2 )    ;
% else
%      Dx(2)=(2/3/m*(-K*(x(1)-l0)*x(2) -Tau* DTheta)       -  x(1)^2*DTheta*D2Theta - x(1)*x(2)*DTheta^2)    /x(2);
% Dx(2)=(2/1/m*(-K*(x(1)-l0)*x(2) + Tau* DTheta)       -  x(1)^2*DTheta*D2Theta - x(1)*x(2)*DTheta^2)    /x(2);
% end


% Dx(3)= Tau/K/(x(3)-l0);

end