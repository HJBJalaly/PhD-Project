function [sys,x0]=Controller(t,x,u, flag,ts)

global PinLock
global t1
global X1pos
global X2pos

x0=[];

switch flag,
  
  

  % --------------
  % Initialization
  % --------------
  case 0,
      
      PinLock=1;
      X1pos=0;
      X2pos=0;
      
      sys=[0 0 3 4 0 0];
    
    % set initial conditions (no conditions)
    x0=[];

    % get actual time
    t1 = now;


      
  % -----------------------------------------------------------------
  % Modification
  % ------------
  case 2, 

    
    % no state return
    sys=[];
    
%     if(PinLock==1)
%         if( u(4) > 0)
%             if(u(3)>u(1))
%                 PinLock=2;
%                 X2pos=u(3);
%             end
%         end
% %     else
%         if( u(2) > 0)
%             if(u(1)>u(3))
%                 PinLock=1;
%                 X1pos=u(1);
%             end
%         end
%     end



  % -------------------------
  % Return Values to Simulink
  % -------------------------
  case 3,     
  
    % no values to return
    sys=[PinLock X1pos X2pos];

 
    
  % ---------------------------
  % Calculate Next Calling Time
  % ---------------------------
  case 4,
    % get number of samples 
    samples = t / ts;
  
    % calculate next calling time
  	sys = (1 + floor( samples + 1e-13*(1+samples) ) )  * ts;



  % -----------------
  % End Of Simulation
  % -----------------
  
  case 9,
    
    % clean up
    sys=[PinLock];
    
    %Zeit erneut messen
    t2 = now;
    fprintf('\nsimulation time: %s\n',datestr(t2-t1, 13));
    


  otherwise
    error(['Unhandled flag = ',num2str(flag)]); % flag error handling


end %switch
