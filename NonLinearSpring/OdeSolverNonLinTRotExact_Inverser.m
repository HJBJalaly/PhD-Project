function Dx=OdeSolverNonLinTRotExact_Inverser(t,x,tau,Theta,L,DL,D2L,Time,K,l0,m)
% Used in NonlinearSpringVerification function
Dx=zeros(2,1);

% Tau   =interp1(Theta,tau,x(1));
Tau   =interp1(Time,tau,t);
ll    =interp1(Theta,L  ,x(1));
dll   =interp1(Theta,DL ,x(1));
d2ll  =interp1(Theta,D2L,x(1));


Dx(1)=x(2);
Dx(2)=(    2/3/m*(  -K*(ll-l0)*dll + Tau* x(2)  )     -  dll*d2ll - ll*dll*x(2)^2 )    /ll^2/x(2);

if(any(isnan(Dx)))
    1;
end

end