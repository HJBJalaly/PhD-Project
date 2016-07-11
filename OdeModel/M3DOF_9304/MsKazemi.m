n=2;
syms a11;
syms a12;
syms a21;
syms a22;
A =([ [a11, a12]; [a21,a22]]);
k=272.2;
T_out=[];
T_sup=[];
p=[];

for j=1:n

    T_out_Temp=input('what is the vector of T_out:');  % T_out_Temp is vertical vector


    T_out(:,end+1)=T_out_Temp;
end
for j1=1:n

    T_sup_Temp=input('what is the vector of T_sup:');  % T_out_Temp is vertical vector

    T_sup(:,end+1)=T_sup_Temp;
end
for j2=1:n

    p_Temp=input('what is the vector of p:');  % T_out_Temp is vertical vector


  p(:,end+1)=p_Temp;
end

%%
xx=solve(k*(T_out-T_sup)-A*k*(T_out- T_sup) - p);
yy=double(struct2array( xx))

       
V=k*(T_out-T_sup);
XX=p/V;
AA=eye(2,2)-XX