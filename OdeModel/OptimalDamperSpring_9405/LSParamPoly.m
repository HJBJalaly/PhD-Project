function  [Beta,Theta] = LSParamPoly(Qq,Dq,Torque,rP,rD,Landa)

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


if(nargin==6)

    Q=[];
    D2Q=[];

    for i=1:length(Qq)
        Q(end+1,:) = Qq(i).^(rP:-1:0)';
        D2Q(end+1,:) = ([Qq(i).^(rP-2:-1:0) 0 0].*(rP:-1:0).*(rP-1:-1:-1))';
    end
    

    % 
    
    % for i=1 to n-2
    A=(Q'*Q+Landa*D2Q'*D2Q);

    Y=Q'*Torque;
    
    Beta=A\Y;

end


end