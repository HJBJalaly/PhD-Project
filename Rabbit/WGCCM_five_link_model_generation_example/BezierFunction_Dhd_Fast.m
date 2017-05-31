function DH=BezierFunction_Dhd_Fast(Ma,Theta,Alfa,Theta_plus,Theta_minus)

s=(Theta-Theta_plus)/(Theta_minus-Theta_plus);
DH=zeros(4,1);
% for i=1:4 % number of output
    
    k=0;
    DH=DH+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
               -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
     
    k=1;
    DH=DH +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                      k*s.^(k-1).*(1-s).^(Ma-k) -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
          
    for k=2:Ma-2
        DH=DH+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
              k*s.^(k-1).*(1-s).^(Ma-k) -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
    end
    
    k=Ma-1;
    DH=DH+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                   k*s.^(k-1).*(1-s).^(Ma-k) -(Ma-k)*s.^(k).*(1-s).^(Ma-k-1)  );
          
    k=Ma;
    DH=DH +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*( ...
                   k*s.^(k-1).*(1-s).^(Ma-k)   );
end
        