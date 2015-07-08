function Dl=OdeSolverNonLinTRot(theta,l,tau,ThetaRange,K,l0)

Tau=interp1(ThetaRange,tau,theta);

Dl= Tau/K/(l-l0);
end