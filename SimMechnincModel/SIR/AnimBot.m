function [sys,x0]=AnimBot(t,x,u, flag,ts); 

% AnimBot.m:      S-function for the animation of the FumJuBot Model 
%																						                          
%																					                             
% Last modified:   September, 4th, 2004                               
%            by:   H. Geyer  				                                  
%





% **************** %
% DECLARATION PART %
% **************** %

% global variables
global   FigHndl         ...        % figure handle
         t1                         % time
  
% view window 
ViewWin = 3;  %[m]

% figure name (identifier)
FigName= 'SIR leg Simulation Model';





% ************ %
% PROGRAM PART %
% ************ %

switch flag,
  
  

  % --------------
  % Initialization
  % --------------
  
  case 0,
     
    
    % Initialize Figure
    % -----------------
    
    % initialize animation figure
    AnimBot_Init(FigName);
    
    % store figure handle for repeated access
    FigHndl = findobj('Type', 'figure',  'Name', FigName);

    % set figure axis range 
    axis([-ViewWin/2 ViewWin/2 -2 1]);
    
    
    % Initialize Plot Handle, i.e. create plot dummy
    % ----------------------------------------------
    
    
    % Annotation: the Simulink inputs u(i) are not present at
    % flag = 0. Hence, the plot dummy must be initiated with
    % all values set to an arbitrary value (zero), i.e. a
    % dummy is created.
    d = [0; 0];
        
    MOHndl = plot(d, d, 'Color', [0.4 0.4 1], 'EraseMode','xor', 'LineWidth', 12);
    GRHndl = plot(d, d, 'Color', 0.3*[1 1 1], 'EraseMode','xor', 'LineWidth', 2);
    % pin1 
    P1Hndl = plot(d, d, 'r', 'EraseMode','xor', ...
                    'MarkerFaceColor',[0 0 1 ],'MarkerEdgeColor',[0 0 1],...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');
    % Spring 1
    S1Hndl = plot(d, d, 'g', 'EraseMode','xor', 'LineWidth',2);
    % Link 1
    L1Hndl = plot(d, d, 'y', 'EraseMode','xor', 'LineWidth',4);
    L2Hndl = plot(d, d, 'y', 'EraseMode','xor', 'LineWidth',4);
    S2Hndl = plot(d, d, 'g', 'EraseMode','xor', 'LineWidth',2);
    M4Hndl = plot(d, d, 'r', 'EraseMode','xor', ...
                    'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1],...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');
    
    % link handles to current figure and to each other
    set(gca,   'UserData', P1Hndl);
    set(P1Hndl,'UserData', S1Hndl);
    set(S1Hndl,'UserData', L1Hndl);
    set(L1Hndl,'UserData', L2Hndl);
    set(L2Hndl,'UserData', S2Hndl);
    set(S2Hndl,'UserData', M4Hndl);
    set(M4Hndl,'UserData', GRHndl);
    set(GRHndl,'UserData', MOHndl);
   
    % set IO-data: .  .  .  number of Simulink "u(i)" - inputs  .  .
    sys=[0 0 0 12 0 0];
    
    % set initial conditions (no conditions)
    x0=[];

    % get actual time
    t1 = now;


  % ------------
  % Modification
  % ------------
  
  case 2, 

    % search root for FigHndl
    if any( get(0,'Children') == FigHndl )
      
      % check handle validity 
      if strcmp(  get( FigHndl,'Name' ), FigName  )
        
        % set actual figure to handle
        set(0, 'CurrentFigure', FigHndl);
        
        
        % Check whether model in view window 
        % and readjust if not
        % ----------------------------------
        
        % get axis limits
        XLimits = get(gca, 'XLim');
        
        % get min and max of 
        minX = min( u(1:2:end-1) );
        maxX = max( u(1:2:end-1) );
        
        % check right border
        if XLimits(2) < ( minX + ViewWin*1/10 )
          
          set(gca, 'XLim', [minX - ViewWin*1/10  minX + ViewWin*9/10]);
        
        % check right border
        elseif XLimits(1) > ( maxX - ViewWin/10 )
          
          set(gca, 'XLim', [maxX - ViewWin*9/10  maxX + ViewWin*1/10]);
        
        end
         

        % Refresh Plot
        % ------------
        
        % get plot handles
        P1Hndl = get(gca,    'UserData');
        S1Hndl = get(P1Hndl, 'UserData');
        L1Hndl = get(S1Hndl, 'UserData');
        L2Hndl = get(L1Hndl, 'UserData');
        S2Hndl = get(L2Hndl, 'UserData');
        M4Hndl = get(S2Hndl, 'UserData');
        GRHndl = get(M4Hndl, 'UserData');
        MOHndl = get(GRHndl, 'UserData');

        % set new plot values     
        set(GRHndl,'XData', [XLimits(1) XLimits(2)], 'YData', [0 0]);
        set(P1Hndl,'XData', u(1), 'YData', u(2));
        set(S1Hndl,'XData', [u(1) u( 3)], 'YData', [u(2) u( 4)]);
        set(L1Hndl,'XData', [u( 3) u( 5)], 'YData', [u( 4) u( 6)]);
        set(L2Hndl,'XData', [u( 7) u(9)], 'YData', [u(8) u(10)]);
        set(S2Hndl,'XData', [u(9) u(11)], 'YData', [u(10) u(12)]);
        set(M4Hndl,'XData', u(11), 'YData', u(12));

      end 
    end

    % no state return
    sys = []; 



  % -------------------------
  % Return Values to Simulink
  % -------------------------
  
  case 3,     
  
    % no values to return
    sys = []; 

 
    
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
    sys = []; 
    
    %Zeit erneut messen
    t2 = now;
    fprintf('\n simulation time: %s \n\a',datestr(t2-t1, 13));
    


  otherwise
    error(['Unhandled flag = ',num2str(flag)]); % flag error handling


end %switch


% annotation(FigHndl,'textbox',...
%     [0.1 0.7 0.1 0.01],...
%     'String',{num2str(t)},...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 0 0],...
%     'Color',[1 0 0]);


% ************* %
% FUNCTION PART %
% ************* %

function AnimBot_Init(namestr)

% -----------------
% Initialize Figure
% -----------------

% check whether figure exists already
[existFlag,figNumber] = figflag(namestr);

% if not, initialize figure
if ~existFlag,
   
  % define figure element
  h0 = figure( ...
       'Tag',                          namestr, ...
       'Name',                         namestr, ...
       'NumberTitle',                    'off', ...
       'BackingStore',                   'off', ...
       'MenuBar',                       'none', ...
		'Color',                        [0 0 0], ...
       'Position',     [100   300   900   300]);
     
     
     
  % define axes element
  h1 = axes( ...
       'Parent',                            h0, ...
       'Tag',                           'axes', ...    
       'Units',                   'normalized', ...
       'Position',       [0.05 0.05 0.92 0.92], ...
       'FontSize',                          8);
       
end %if ~existflag



% ----------------------------------
% Reset Figure to Simulation Default
% ----------------------------------

% reset axes to default properties
cla reset;

% change some properties
set(gca, 'DrawMode',   'fast', ...
         'Visible',      'on', ...
         'Color',     [0 0 0], ...
         'XColor',    [1 1 1], ...
				 'YColor',    [1 1 1]);

axis on;
axis image;

hold on;


