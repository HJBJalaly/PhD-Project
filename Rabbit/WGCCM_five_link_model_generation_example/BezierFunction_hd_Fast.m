function H=BezierFunction_hd_Fast(Ma,Theta,Alfa,Theta_plus,Theta_minus)

s=(Theta-Theta_plus)/(Theta_minus-Theta_plus);
H=zeros(4,1);
    
    k=0;
    H=H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
     
    k=1;
    H=H   +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
          
    for k=2:Ma-2
        H=H+Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
    end
    
    k=Ma-1;
    H=H   +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
          
    k=Ma;
    H=H   +Alfa(:,k+1)*factorial(Ma)/factorial(k)/factorial(Ma-k)*s.^k.*(1-s).^(Ma-k);
end
        