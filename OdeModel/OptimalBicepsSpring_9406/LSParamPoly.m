function  [Beta,Theta] = LSParamPoly(Qq,Torque,rP)

QQ=[];

for i=1:length(Qq)
     QQ(i,:) = Qq(i).^(rP:-1:0)';
end


Beta=((QQ'*QQ))\( QQ'*Torque);

end