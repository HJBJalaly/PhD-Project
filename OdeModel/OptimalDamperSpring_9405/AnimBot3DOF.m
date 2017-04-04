function AnimBot3DOF(T,Y,L)

FigName= '3DoF Manipulator Simulation Model';
FramesFolder = './ImageExmaple';

% Initialize Figure
% -----------------

% initialize animation figure
AnimBot_Init(FigName);

% store figure handle for repeated access
FigHndl = findobj('Type', 'figure',  'Name', FigName);

% set figure axis range
ViewWin = 2*L;  %[m]
axis([-ViewWin ViewWin -.5*L ViewWin]);


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
% Link 1
L1Hndl = plot(d, d, 'y', 'EraseMode','normal', 'LineWidth',4);
% pin2
P2Hndl = plot(d, d, 'r', 'EraseMode','normal', ...
    'MarkerFaceColor',[0 1 1 ],'MarkerEdgeColor',[0 1 1],...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
% Link 2
L2Hndl = plot(d, d, 'y', 'EraseMode','normal', 'LineWidth',4);

% EF
EFHndl = plot(d, d, 'r', 'EraseMode','normal', ...
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
set(P1Hndl,'UserData', L1Hndl);
set(L1Hndl,'UserData', P2Hndl);
set(P2Hndl,'UserData', L2Hndl);
set(L2Hndl,'UserData', EFHndl);
set(EFHndl,'UserData', GRHndl);
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
        
        
%         % ------------
%         
%         % get plot handles
%         P1Hndl = get(gca,    'UserData');
%         L1Hndl = get(P1Hndl, 'UserData');
%         P2Hndl = get(L1Hndl, 'UserData');
%         L2Hndl = get(P2Hndl, 'UserData');
%         EFHndl = get(L2Hndl, 'UserData');
%         GRHndl = get(EFHndl, 'UserData');
%         MOHndl = get(GRHndl, 'UserData');
%         TrackHandl = get(MOHndl, 'UserData');
%         
        % set new plot values

        Pin1Pos=[0 0];
        Tracker=[];
%         frame=0;
        for i=1:1:length(T)-1
%             frame=frame+1;
            Link1PosA =[Pin1Pos(1)                   , Pin1Pos(2)];
            Link1PosB =[Pin1Pos(1)+(L)*cos(Y(i,1))   , Pin1Pos(2)+(L)*sin(Y(i,1))];
            
            Pin2Pos   = Link1PosB;
            Link2PosA = Pin2Pos;
            Link2PosB =[Pin2Pos(1)+(L)*cos(Y(i,2)+Y(i,1))  , Pin2Pos(2)+(L)*sin(Y(i,2)+Y(i,1))];
            
            
            EFPos     =Link2PosB;
            
            Tracker=[Tracker;EFPos];

            set(GRHndl,'XData', [XLimits(1) XLimits(2)]      , 'YData', [0 0]);
            set(P1Hndl,'XData',  Pin1Pos(1)                  , 'YData', Pin1Pos(2));
            set(L1Hndl,'XData', [Link1PosA(1) Link1PosB(1)]  , 'YData', [Link1PosA(2) Link1PosB(2)]);
            set(P2Hndl,'XData',  Pin2Pos(1)                  , 'YData', Pin2Pos(2));
            set(L2Hndl,'XData', [Link2PosA(1) Link2PosB(1)]  , 'YData', [Link2PosA(2) Link2PosB(2)]);

            set(EFHndl,'XData', EFPos(1)                     , 'YData', EFPos(2));
            
            set(TrackHandl,'XData', Tracker(:,1), 'YData', Tracker(:,2));

            drawnow;
%             SavePdfFast(sprintf('%s/Circle_%02d', FramesFolder,frame))            
%             pause((T(i+1)-T(i)))
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