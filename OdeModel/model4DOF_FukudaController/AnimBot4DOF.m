function AnimBot4DOF(T,Y,Pin1Pos,L,d0)

FigName= 'SIR leg Simulation Model';


% Initialize Figure
% -----------------

% initialize animation figure
AnimBot_Init(FigName);

% store figure handle for repeated access
FigHndl = findobj('Type', 'figure',  'Name', FigName);

% set figure axis range
ViewWin = 3;  %[m]
axis([-ViewWin/2 ViewWin/2 -2 1]);


% Initialize Plot Handle, i.e. create plot dummy
% ----------------------------------------------


% Annotation: the Simulink inputs u(i) are not present at
% flag = 0. Hence, the plot dummy must be initiated with
% all values set to an arbitrary value (zero), i.e. a
% dummy is created.
d = [0; 0];

MOHndl = plot(d, d, 'Color', [0.4 0.4 1], 'EraseMode','normal', 'LineWidth', 12);
GRHndl = plot(d, d, 'Color', 0.3*[1 1 1], 'EraseMode','normal', 'LineWidth', 2);
% pin1
P1Hndl = plot(d, d, 'r', 'EraseMode','normal', ...
    'MarkerFaceColor',[0 0 1 ],'MarkerEdgeColor',[0 0 1],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
% Spring 1
S1Hndl = plot(d, d, 'g', 'EraseMode','normal', 'LineWidth',2);
% Link 1
L1Hndl = plot(d, d, 'y', 'EraseMode','normal', 'LineWidth',4);
P2Hndl = plot(d, d, 'r', 'EraseMode','normal', ...
    'MarkerFaceColor',[0 1 1 ],'MarkerEdgeColor',[0 1 1],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
L2Hndl = plot(d, d, 'y', 'EraseMode','normal', 'LineWidth',4);
S2Hndl = plot(d, d, 'g', 'EraseMode','normal', 'LineWidth',2);
M4Hndl = plot(d, d, 'r', 'EraseMode','normal', ...
    'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
TrackHandl= plot(d, d, 'Color',[1 1 1], 'EraseMode','background', ...
    'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1],...
    'MarkerSize',1,...
    'Marker','none',...
    'LineStyle','-.');

% link handles to current figure and to each other
set(gca,   'UserData', P1Hndl);
set(P1Hndl,'UserData', S1Hndl);
set(S1Hndl,'UserData', L1Hndl);
set(L1Hndl,'UserData', P2Hndl);
set(P2Hndl,'UserData', L2Hndl);
set(L2Hndl,'UserData', S2Hndl);
set(S2Hndl,'UserData', M4Hndl);
set(M4Hndl,'UserData', GRHndl);
set(GRHndl,'UserData', MOHndl);
set(MOHndl,'UserData', TrackHandl);


% ------------
% Modification
% ------------


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
%         minX = min( u(1:2:end-1) );
%         maxX = max( u(1:2:end-1) );
%         
%         % check right border
%         if XLimits(2) < ( minX + ViewWin*1/10 )
%             
%             set(gca, 'XLim', [minX - ViewWin*1/10  minX + ViewWin*9/10]);
%             
%             % check right border
%         elseif XLimits(1) > ( maxX - ViewWin/10 )
%             
%             set(gca, 'XLim', [maxX - ViewWin*9/10  maxX + ViewWin*1/10]);
%             
%         end
%         
        
        
        % ------------
        
        % get plot handles
        P1Hndl = get(gca,    'UserData');
        S1Hndl = get(P1Hndl, 'UserData');
        L1Hndl = get(S1Hndl, 'UserData');
        P2Hndl = get(L1Hndl, 'UserData');
        L2Hndl = get(P2Hndl, 'UserData');
        S2Hndl = get(L2Hndl, 'UserData');
        M4Hndl = get(S2Hndl, 'UserData');
        GRHndl = get(M4Hndl, 'UserData');
        MOHndl = get(GRHndl, 'UserData');
        TrackHandl = get(MOHndl, 'UserData');
        
        % set new plot values

        Tracker=[];
        for i=1:length(T)-1
            Link1UpPos=[Pin1Pos(1)+(Y(i,1)+d0)*sin(Y(i,2))     , Pin1Pos(2)-(Y(i,1)+d0)*cos(Y(i,2))];
            Link1DnPos=[Pin1Pos(1)+(Y(i,1)+d0+L)*sin(Y(i,2))   , Pin1Pos(2)-(Y(i,1)+d0+L)*cos(Y(i,2))];
            Link2DnPos= Link1DnPos;
            Pin2Pos   = Link1DnPos;
            Link2UpPos=[Link1DnPos(1)+(L)*sin(Y(i,4)-Y(i,2))  , Link1DnPos(2)+(L)*cos(Y(i,4)-Y(i,2))];
            Pin4Pos   =[Link1DnPos(1)+(L+d0+Y(i,3))*sin(Y(i,4)-Y(i,2))  , Link1DnPos(2)+(L+d0+Y(i,3))*cos(Y(i,4)-Y(i,2))];
            Tracker=[Tracker;Pin4Pos];

            
            
            set(GRHndl,'XData', [XLimits(1) XLimits(2)]      , 'YData', [0 0]);
            set(P1Hndl,'XData',  Pin1Pos(1)                  , 'YData', Pin1Pos(2));
            set(P2Hndl,'XData', Pin2Pos(1)                   , 'YData', Pin2Pos(2));
            set(M4Hndl,'XData', Pin4Pos(1)                   , 'YData', Pin4Pos(2));
            
            set(S1Hndl,'XData', [Pin1Pos(1) Link1UpPos(1)]   , 'YData', [Pin1Pos(2) Link1UpPos(2)]);
            set(L1Hndl,'XData', [Link1UpPos(1) Link1DnPos(1)], 'YData', [Link1UpPos(2) Link1DnPos(2)]);
            set(L2Hndl,'XData', [Link2DnPos(1) Link2UpPos(1)], 'YData', [Link2DnPos(2) Link2UpPos(2)]);
            set(S2Hndl,'XData', [Link2UpPos(1) Pin4Pos(1)]   , 'YData', [Link2UpPos(2) Pin4Pos(2)]);
            
            set(TrackHandl,'XData', Tracker(:,1), 'YData', Tracker(:,2));

            drawnow;
            pause(4*(T(i+1)-T(i)))
        end
    end
end


end

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
end