function  [Beta,Theta] = LSParamPoly(Qq,Dq,Torque,rP,rD)

QQ=[];
DD2=[];
for i=1:length(Qq)
     QQ(i,:) = Qq(i).^(rP:-1:0)';
     DQ(i,:) = Dq(i).^(rD:-1:0)';
end


if(rP>0)
    if(rD>0)
        Temp=(eye(length(DQ))- DQ*((DQ'*DQ)\DQ'));
    else
        Temp=eye(length(DQ));
    end
    Beta=((QQ'*Temp*QQ))\( QQ'*Temp*Torque);
else
    Beta=[0 0];
end


if(rD>0)
    if(rP>0)
        Gamma=(eye(length(DQ))- DQ*((DQ'*DQ)\DQ'));
        Theta=(DQ'*DQ)\( DQ'*Torque-DQ'*QQ*Beta );
    else
        Theta=((DQ'*DQ))\( DQ'*Torque);
    end
else
    Theta=[0 0];
end


end