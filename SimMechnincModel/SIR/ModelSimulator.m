function ModelSimulator()

clear
close all
Initialization();

XInit=0;
YInit=0;
ThInit=-60;  % for Pin 1 (AngPin)
ThetaInc=-50;   % for angle between tw0 link
sFInit1=-0.1;
sFInit2=-0.1;

open('SIR_Temp6')
set_param('SIR_Temp6/L1 teta1  IC3','P1PIC',num2str(XInit));
set_param('SIR_Temp6/L1 teta1  IC3','P2PIC',num2str(YInit));
set_param('SIR_Temp6/L1 teta1  IC3','R1PIC',num2str(ThInit));
set_param('SIR_Temp6/Spring Free lenght  1','P1PIC',num2str(sFInit1));
set_param('SIR_Temp6/Spring Free lenght  2','P1PIC',num2str(sFInit2));
set_param('SIR_Temp6/LinkAngle','R1PIC',num2str(ThetaInc));

sim('SIR_Temp6',[0 10])

for i=1:10
    
    sFInit2=SpringLength1.signals.values(end);
    sFInit1=SpringLength2.signals.values(end);
    
    ThInit=(Theta.signals.values(end)+AngPin1.signals.values(end));  % for Pin 1
    ThetaInc=-Theta.signals.values(end);  % for angle between tw0 link
    
    XInit=M4C.signals.values(end,1)+0.025*sind(abs(ThInit));
    YInit=M4C.signals.values(end,2)+0.01;
    
    set_param('SIR_Temp6/L1 teta1  IC3','P1PIC',num2str(XInit));
    set_param('SIR_Temp6/L1 teta1  IC3','P2PIC',num2str(YInit));
    set_param('SIR_Temp6/L1 teta1  IC3','R1PIC',num2str(ThInit));
    set_param('SIR_Temp6/Spring Free lenght  1','P1PIC',num2str(sFInit1));
    set_param('SIR_Temp6/Spring Free lenght  2','P1PIC',num2str(sFInit2));
    set_param('SIR_Temp6/LinkAngle','R1PIC',num2str(ThetaInc));

    sim('SIR_Temp6',[0 10])
end