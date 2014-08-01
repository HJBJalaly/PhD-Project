g=9.81;

mPin1=0.5;
rPin1=0.025;
sK1=1000;
sFl1=-0.10; % free lenght of spring
mLink1=1;
lenLink1=0.5;

mPin2=0.5;
rPin2=0.025;
sK2=1000;
sFl2=-0.10; % free lenght of spring
mLink2=1;
lenLink2=0.5;


% Ground Interaction
% ------------------

% stiction coefficient
mu_stick = 1.05;

% sliding friction coefficient
mu_slide = .9;

pengain = 400000;
pendamp= 15;

global ReachLine
ReachLine=0;