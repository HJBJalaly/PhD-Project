function  Ainv=SVDBlockInvertor(A,numberB,sizeB,MinSinValue)

% A is a block dignal matrix, which has "numberB" of blocks.
% Each block is a "sizeB" squer matrtix.

Ainv=zeros(numberB*sizeB);

for i=1:numberB
   Ablock=A((i-1)*sizeB+1:(i)*sizeB,(i-1)*sizeB+1:(i)*sizeB);
   [U,S,V]=svd(Ablock);
   [rr,cc]=find(S>MinSinValue);
   Sp=zeros(sizeB);
   Sr=S(1:max(rr),1:max(rr));
   Sp(1:max(rr),1:max(rr))=Sr^-1;
   AblockInv=V*Sp*U';
   Ainv((i-1)*sizeB+1:(i)*sizeB,(i-1)*sizeB+1:(i)*sizeB)=AblockInv;
end




