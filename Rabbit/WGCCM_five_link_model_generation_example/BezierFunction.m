function [H,DH,D2H]=BezierFunction(Theta,Alfa,Theta_plus,Theta_minus)

s=(Theta-Theta_plus)/(Theta_minus-Theta_plus);
Ma=3;
H=zeros(4,1);
DH=zeros(4,1);
D2H=zeros(4,1);
% for i=1:4 % number of output
    
    k=0;
    H=H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
    DH=DH+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
               -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
    D2H=D2H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                  +(Ma-k)*(Ma-k-1)*s.^(k).*(1-s).^(Ma-k-2)  );
     
    k=1;
    H=H   +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
    DH=DH +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                      k*s.^(k-1).*(1-s).^(Ma-k) -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
    D2H=D2H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                     -2*(Ma-k)*k*s.^(k-1).*(1-s).^(Ma-k-1) +(Ma-k)*(Ma-k-1)*s.^(k).*(1-s).^(Ma-k-2)  );
          
    for k=2:Ma-2
        H=H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
        DH=DH+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
              k*s.^(k-1).*(1-s).^(Ma-k) -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
        D2H=D2H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
              k*(k-1)*s.^(k-2).*(1-s).^(Ma-k) -2*(Ma-k)*k*s.^(k-1).*(1-s).^(Ma-k-1) +(Ma-k)*(Ma-k-1)*s.^(k).*(1-s).^(Ma-k-2)  );
    end
    
    k=Ma-1;
    H=H   +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
    DH=DH+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                   k*s.^(k-1).*(1-s).^(Ma-k) -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
    D2H=D2H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                   k*(k-1)*s.^(k-2).*(1-s).^(Ma-k) -2*(Ma-k)*k*s.^(k-1).*(1-s).^(Ma-k-1));
          
    k=Ma;
    H=H   +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
    DH=DH +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                   k*s.^(k-1).*(1-s).^(Ma-k)   );
    D2H=D2H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                   k*(k-1)*s.^(k-2).*(1-s).^(Ma-k)  );
end

% h=[];
% for i=1:4
%     s=0:.02:1;
%     h(i,length(s))=0;
%     for k=0:Ma
%         h(i,:)=h(i,:)+Alfa(i,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
%     end
%     subplot(2,2,i)
%     plot(s,h(i,:))
% end
% %     
        