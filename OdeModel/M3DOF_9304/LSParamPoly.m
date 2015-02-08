function  CoefBLS = LSParamPoly(qtrajectory,TorqueProfile,rU,alpha)

QQ=[];
DD=[];
for i=1:length(qtrajectory)
     QQ(i,:) = qtrajectory(i).^(rU:-1:0)';
     DD(i,:) = ([qtrajectory(i).^(rU-1:-1:0) 0].*(rU:-1:0))';
end


% CoefBLS=((QQ'*QQ)+alpha*(DD'*DD))^-1 * QQ'*TorqueProfile;
CoefBLS=((QQ'*QQ)+alpha*(DD'*DD))\( QQ'*TorqueProfile);

% syms s real;
% 
% Psi=zeros(rU+1,rU+1);
% 
% c_hat =(max(qtrajectory) - min(qtrajectory));
% d_hat = min(qtrajectory);
% 
% for kk=1:rU+1
%     Psi(kk,:)  = [zeros(1,kk-1),sym2poly( (c_hat*s+ d_hat)^(rU+1-(kk))  ) ];
% end
% 
% qtrajectory_hat= ( qtrajectory-d_hat)/c_hat/2+.2;
% 
% for i=1:length(qtrajectory)
%      QQ_hat(i,:) = qtrajectory_hat(i).^(rU:-1:0)';
%      DD_hat(i,:) = ([qtrajectory_hat(i).^(rU-1:-1:0) 0].*(rU:-1:0))';
% end
% 
% CoefBLS_hat=((QQ_hat'*QQ_hat)+alpha*(DD_hat'*DD_hat))\( QQ_hat'*TorqueProfile);

end