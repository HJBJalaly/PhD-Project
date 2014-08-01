function Ocilator()
home
[T,X]=ode23s(@(t,x)System(t,x),[0:.01:100],[3.5 0]);
subplot(1,2,1)
plot(T,X(:,1))
subplot(1,2,2)
plot(T,0.21*(X(:,1)-0.75*pi).*(X(:,1)-0.25*pi).*(X(:,1)-1.25*pi)+2.5)
end



function Dx=System(t,x)

f=0.21*(x(1)-0.75*pi)*(x(1)-0.25*pi)*(x(1)-1.25*pi)+2.5;

Dx=zeros(2,1);
Dx(1)=x(2);
Dx(2)=-f;
end

