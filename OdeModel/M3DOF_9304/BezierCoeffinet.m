function  [CoefP,CoefB]= BezierCoeffinet(Tspan,Sspan,Degree)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Bezier=0;
Alfha=sym('Alp',[Degree+1,1]);
for i=0:Degree
   Bezier=Bezier + Alfha(i+1)*factorial(Degree)/(factorial(i)*factorial(Degree-i))*(-1)^(Degree-i)*poly([ones(1,Degree-i) zeros(1,i)]);
end


TspanNor=(Tspan-min(Tspan))/(max(Tspan)-min(Tspan));
[CoefP,Po] = polyfit(Tspan,Sspan,Degree);
[CoefPB,Po] = polyfit(TspanNor,Sspan,Degree);
CoefB=double(struct2array(solve(Bezier-CoefPB)));

end