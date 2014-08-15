function [rr_rest,theta_rest] = InitCondinRest(Space,LL,d0,mL,mp,Ks,gg)

S=Space/2;
L=LL+d0;
M=(mp/2+mL)*gg/Ks;

syms x
x=solve(S*cot(x)-L*cos(x)-M);

theta_rest=rad2deg(double(abs(x)));
rr_rest=S/sind(theta_rest) - L;

end

