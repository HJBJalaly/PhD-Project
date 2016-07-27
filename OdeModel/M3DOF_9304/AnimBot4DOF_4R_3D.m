function AnimBot4DOF_4R_3D(T,Y,LL1,LL2,LL3,LL4)

FigName= '4DoF Manipulator Simulation Model';


% Initialize Figure
% -----------------

% initialize animation figure
AnimBot_Init(FigName);

% store figure handle for repeated access
FigHndl = findobj('Type', 'figure',  'Name', FigName);

% set figure axis range
ViewWin = 1;  %[m]
axis([-ViewWin ViewWin -ViewWin ViewWin -ViewWin ViewWin ]);


% Initialize Plot Handle, i.e. create plot dummy
% ----------------------------------------------


% Annotation: the Simulink inputs u(i) are not present at
% flag = 0. Hence, the plot dummy must be initiated with
% all values set to an arbitrary value (zero), i.e. a
% dummy is created.
d = [0; 0; 0];

SUB={1 2,3 ,4};
for sub=1:4;
    
    MOHndl{sub} = plot3(d, d,d, 'Color', [0.4 0.4 1], 'EraseMode','normal', 'LineWidth', 12);
    GRHndl{sub} = plot3(d, d,d, 'Color', 0.3*[1 1 1], 'EraseMode','normal', 'LineWidth', 2);
    % pin1
    P1Hndl{sub} = plot3(d, d,d, 'r', 'EraseMode','normal', ...
        'MarkerFaceColor',[0 0 1 ],'MarkerEdgeColor',[0 0 1],...
        'MarkerSize',10,...
        'Marker','o',...
        'LineStyle','none');
    % Link 1
    L1Hndl{sub} = plot3(d, d,d, 'y', 'EraseMode','normal', 'LineWidth',4);
    % pin2
    P2Hndl{sub} = plot3(d, d,d, 'r', 'EraseMode','normal', ...
        'MarkerFaceColor',[0 1 1 ],'MarkerEdgeColor',[0 1 1],...
        'MarkerSize',10,...
        'Marker','o',...
        'LineStyle','none');
    % Link 2
    L2Hndl{sub} = plot3(d, d,d, 'y', 'EraseMode','normal', 'LineWidth',4);
    % pin3
    P3Hndl{sub} = plot3(d, d,d, 'r', 'EraseMode','normal', ...
        'MarkerFaceColor',[0 1 1 ],'MarkerEdgeColor',[0 1 1],...
        'MarkerSize',10,...
        'Marker','o',...
        'LineStyle','none');
    % Link 3
    L3Hndl{sub} = plot3(d, d,d, 'y', 'EraseMode','normal', 'LineWidth',4);

    % pin4
    P4Hndl{sub} = plot3(d, d,d, 'r', 'EraseMode','normal', ...
        'MarkerFaceColor',[0 1 1 ],'MarkerEdgeColor',[0 1 1],...
        'MarkerSize',10,...
        'Marker','o',...
        'LineStyle','none');
    % Link 4
    L4Hndl{sub} = plot3(d, d,d, 'y', 'EraseMode','normal', 'LineWidth',4);

    % EF
    EFHndl{sub} = plot3(d, d,d, 'r', 'EraseMode','normal', ...
        'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1],...
        'MarkerSize',10,...
        'Marker','o',...
        'LineStyle','none');
    TrackHandl{sub}= plot3(d, d,d, 'Color',[1 1 1], 'EraseMode','background', ...
        'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1],...
        'MarkerSize',1,...
        'Marker','none',...
        'LineStyle','-.');

    % link handles to current figure and to each other
    subplot(2,2,SUB{sub})
    set(gca,   'UserData', P1Hndl{sub});
    set(P1Hndl{sub},'UserData', L1Hndl{sub});
    set(L1Hndl{sub},'UserData', P2Hndl{sub});
    set(P2Hndl{sub},'UserData', L2Hndl{sub});
    set(L2Hndl{sub},'UserData', P3Hndl{sub});
    set(P3Hndl{sub},'UserData', L3Hndl{sub});
    set(L3Hndl{sub},'UserData', P4Hndl{sub});
    set(P4Hndl{sub},'UserData', L4Hndl{sub});
    set(L4Hndl{sub},'UserData', EFHndl{sub});
    set(EFHndl{sub},'UserData', GRHndl{sub});
    set(GRHndl{sub},'UserData', MOHndl{sub});
    set(MOHndl{sub},'UserData', TrackHandl{sub});

end
% ------------
% Modification
% ------------


% search root for FigHndl
if any( get(0,'Children') == FigHndl )
    
    % check handle validity
    if strcmp(  get( FigHndl,'Name' ), FigName  )
        
        % set actual figure to handle
        set(0, 'CurrentFigure', FigHndl);
        
                
        % get axis limits
        XLimits = get(gca, 'XLim');
        
        % set new plot values

        Pin1Pos=[0 0 0];
        Tracker=[];
        for i=1:4:length(T)-1
            q1=Y(i,1);
            q2=Y(i,2);
            q3=Y(i,3);
            q4=Y(i,4);
            
            Link1PosA =Pin1Pos;
            Link1PosB =[ 0, 0, LL1];
            
            Pin2Pos   = Link1PosB;
            Link2PosA = Pin2Pos;
            Link2PosB = Pin2Pos+[ LL2*cos(q1)*sin(q2), LL2*sin(q1)*sin(q2), LL2*cos(q2)];
            
            Pin3Pos   = Link2PosB;
            Link3PosA = Pin3Pos;
            Link3PosB = Pin3Pos +[ LL3*sin(q2 + q3)*cos(q1), LL3*sin(q2 + q3)*sin(q1), LL3*cos(q2 + q3)];
            
            Pin4Pos   = Link3PosB;
            Link4PosA = Pin4Pos;
            Link4PosB = Pin4Pos +[ LL4*sin(q2 + q3 + q4)*cos(q1), LL4*sin(q2 + q3 + q4)*sin(q1), LL4*cos(q2 + q3 + q4)];
            
            EFPos     =Link4PosB;
            
            Tracker=[Tracker;EFPos];

            for sub=1:4
                set(0, 'CurrentFigure', FigHndl);
                subplot(2,2,SUB{sub})
                set(GRHndl{sub},'XData', [XLimits(1) XLimits(2) ]     , 'YData', [0 0]                      , 'ZData', [0 0]);
                set(P1Hndl{sub},'XData',  Pin1Pos(1)                  , 'YData', Pin1Pos(2)                 , 'ZData', Pin1Pos(3));
                set(L1Hndl{sub},'XData', [Link1PosA(1) Link1PosB(1)]  , 'YData', [Link1PosA(2) Link1PosB(2)], 'ZData', [Link1PosA(3) Link1PosB(3)]);
                set(P2Hndl{sub},'XData',  Pin2Pos(1)                  , 'YData'    , Pin2Pos(2)             , 'ZData', Pin2Pos(3));
                set(L2Hndl{sub},'XData', [Link2PosA(1) Link2PosB(1)]  , 'YData', [Link2PosA(2) Link2PosB(2)], 'ZData', [Link2PosA(3) Link2PosB(3)]);
                set(P3Hndl{sub},'XData',  Pin3Pos(1)                  , 'YData'    , Pin3Pos(2)             , 'ZData', Pin3Pos(3));
                set(L3Hndl{sub},'XData', [Link3PosA(1) Link3PosB(1)]  , 'YData', [Link3PosA(2) Link3PosB(2)], 'ZData', [Link3PosA(3) Link3PosB(3)]);
                set(P4Hndl{sub},'XData',  Pin4Pos(1)                  , 'YData'    , Pin4Pos(2)             , 'ZData', Pin4Pos(3));
                set(L4Hndl{sub},'XData', [Link4PosA(1) Link4PosB(1)]  , 'YData', [Link4PosA(2) Link4PosB(2)], 'ZData', [Link4PosA(3) Link4PosB(3)]);

                set(EFHndl{sub},'XData', EFPos(1)        , 'YData', EFPos(2)    , 'ZData', EFPos(3));
                set(TrackHandl{sub},'XData', Tracker(:,1), 'YData', Tracker(:,2), 'ZData', Tracker(:,3));

                drawnow;
            end
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
if existFlag,
    clf

    figure( figNumber)
    set(gcf,...
        'NumberTitle',                    'off', ...
        'BackingStore',                   'off', ...
        'MenuBar',                       'none', ...
        'Color',                        [0 0 0], ...
        'Position',     [100   300   900   600]);
else
    figure( ...
        'Tag',                          namestr, ...
        'Name',                         namestr, ...
        'NumberTitle',                    'off', ...
        'BackingStore',                   'off', ...
        'MenuBar',                       'none', ...
        'Color',                        [0 0 0], ...
        'Position',     [100   300   900   600]);

end

% change some properties
subplot(2,2,1)    
plot3(0,0,0);
set(gca, 'DrawMode',   'fast', ...
    'Visible',      'on', ...
    'Color',     [0 0 0], ...
    'XColor',    [1 1 1], ...
    'YColor',    [1 1 1], ...
    'ZColor',    [1 1 1]);
axis on;
axis equal;
hold on;

subplot(2,2,2)    
plot3(0,0,0);
xlabel('x')
ylabel('y')
set(gca,'CameraPosition', [0 0 17.3205]) % x y
set(gca, 'DrawMode',   'fast', ...
    'Visible',      'on', ...
    'Color',     [0 0 0], ...
    'XColor',    [1 1 1], ...
    'YColor',    [1 1 1], ...
    'ZColor',    [1 1 1]);
axis on;
axis equal;
hold on;


subplot(2,2,3)    
plot3(0,0,0);
xlabel('x')
zlabel('z')
set(gca,'CameraPosition', [0 -17.3205 0]) % x z
set(gca, 'DrawMode',   'fast', ...
    'Visible',      'on', ...
    'Color',     [0 0 0], ...
    'XColor',    [1 1 1], ...
    'YColor',    [1 1 1], ...
    'ZColor',    [1 1 1]);
axis on;
axis equal;
hold on;

subplot(2,2,4)    
plot3(0,0,0);
set(gca,'CameraPosition', [17.3205 0 0])  % y z
ylabel('y')
zlabel('z')
set(gca, 'DrawMode',   'fast', ...
    'Visible',      'on', ...
    'Color',     [0 0 0], ...
    'XColor',    [1 1 1], ...
    'YColor',    [1 1 1], ...
    'ZColor',    [1 1 1]);
axis on;
axis equal;
hold on;
end
