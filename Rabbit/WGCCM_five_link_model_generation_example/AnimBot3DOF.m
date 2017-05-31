function AnimBot3DOF(time,q,T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,PaTi)

FigName= 'Rabbit Simulation Model';
% FramesFolder = './ImageExmaple';


% Initialize Figure
% -----------------

% initialize animation figure
AnimBot_Init(FigName);

% store figure handle for repeated access
FigHndl = findobj('Type', 'figure',  'Name', FigName);

% set figure axis range
ViewWin = 8*Lc_fem;  %[m]
axis([-1 10 -ViewWin ViewWin]);


% Initialize Plot Handle, i.e. create plot dummy
% ----------------------------------------------


% Annotation: the Simulink inputs u(i) are not present at
% flag = 0. Hence, the plot dummy must be initiated with
% all values set to an arbitrary value (zero), i.e. a
% dummy is created.
d = [0; 0];

MOHndl = plot(d, d, 'Color', [0.4 0.4 1], 'EraseMode','normal', 'LineWidth', 12);
GRHndl = plot(d, d, 'Color', 0.3*[1 1 1], 'EraseMode','normal', 'LineWidth', 2);
% foot 1
foot1Hndl = plot(d, d, 'b', 'EraseMode','normal', ...
                    'MarkerFaceColor','b','MarkerEdgeColor','b',...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');
% tibia 1
Tibia1Hndl = plot(d, d, 'b', 'EraseMode','normal', 'LineWidth',4);
% knee 1
Knee1Hndl = plot(d, d, 'b', 'EraseMode','normal', ...
                    'MarkerFaceColor','b','MarkerEdgeColor','b',...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');
% femur 1 
Femur1Hndl = plot(d, d, 'b', 'EraseMode','normal', 'LineWidth',4);
% hip 
Hip_Hndl = plot(d, d, 'r', 'EraseMode','normal', ...
                    'MarkerFaceColor','r','MarkerEdgeColor','r',...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');
% Torso
Torso_Hndl = plot(d, d, 'y', 'EraseMode','normal', 'LineWidth',4);
% femur 2
Femur2Hndl = plot(d, d, 'g', 'EraseMode','normal', 'LineWidth',4);

% knee 2
Knee2Hndl = plot(d, d, 'g', 'EraseMode','normal', ...
                    'MarkerFaceColor','g','MarkerEdgeColor','g',...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');

% tibia 2
Tibia2Hndl = plot(d, d, 'g', 'EraseMode','normal', 'LineWidth',4);
% foot 2
foot2Hndl = plot(d, d, 'g', 'EraseMode','normal', ...
                    'MarkerFaceColor','g','MarkerEdgeColor','g',...
                    'MarkerSize',10,...
                    'Marker','o',...
                    'LineStyle','none');
                
                
TrackHandl= plot(d, d, 'Color',[1 1 1], 'EraseMode','background', ...
    'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1],...
    'MarkerSize',1,...
    'Marker','none',...
    'LineStyle','-.');

% link handles to current figure and to each other
set(gca,   'UserData', foot1Hndl);
set(foot1Hndl,'UserData', Tibia1Hndl);
set(Tibia1Hndl,'UserData', Knee1Hndl);
set(Knee1Hndl,'UserData', Femur1Hndl);
set(Femur1Hndl,'UserData', Hip_Hndl);
set(Hip_Hndl,'UserData', Torso_Hndl);
set(Torso_Hndl,'UserData', Femur2Hndl);
set(Femur2Hndl,'UserData', Knee2Hndl);
set(Knee2Hndl,'UserData', Tibia2Hndl);
set(Tibia2Hndl,'UserData', foot2Hndl);
set(foot2Hndl,'UserData', GRHndl);
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
        
        TT = [ 1   0   0   0  -1;
              0   1   0   0  -1;
             -1   0   1   0   0;
              0  -1   0   1   0;
              0   0   0   0   1];

    % old absolute coordinates in terms of relative coords
%     q_new  = inv(TT)*q;
    


    Tracker=[];
    ImpactCounter=1;
    p_foot1    =[0 0]';
        
%     frame=0;
    for i=1:1:length(time)-1
        
        q_fem1  = q(i,1)+q(i,5);
        q_fem2  = q(i,2)+q(i,5);
        q_tib1  = q(i,1)+q(i,3)+q(i,5);
        q_tib2  = q(i,2)+q(i,4)+q(i,5);
        q_torso = q(i,5);

        p_cen_tib1 =p_foot1 + (L_tib-Lc_tib   )*[ sin(pi-q_tib1), cos(pi-q_tib1)]';
        p_knee1    =p_foot1 + (L_tib          )*[ sin(pi-q_tib1), cos(pi-q_tib1)]';
        p_cen_fem1 =p_knee1+ (L_fem-Lc_fem    )*[ sin(pi-q_fem1), cos(pi-q_fem1)]';
        p_hip      =p_knee1+ (L_fem           )*[ sin(pi-q_fem1), cos(pi-q_fem1)]';
        p_cen_torso=p_hip  + (L_torso-Lc_torso)*[-sin(q_torso)  , cos(q_torso)  ]';
        p_head     =p_hip  + (L_torso         )*[-sin(q_torso)  , cos(q_torso)  ]';
        p_cen_fem2 =p_hip  + (Lc_fem          )*[ sin(q_fem2-pi),-cos(q_fem2-pi)]';
        p_knee2    =p_hip  + (L_fem           )*[ sin(q_fem2-pi),-cos(q_fem2-pi)]';
        p_cen_tib2 =p_knee2+ (Lc_tib          )*[ sin(q_tib2-pi),-cos(pi-q_tib2)]';
        p_foot2    =p_knee2+ (L_tib           )*[ sin(q_tib2-pi),-cos(pi-q_tib2)]';

        set(GRHndl,'XData',     [XLimits(1) XLimits(2)]     , 'YData', [0 0]);
        set(foot1Hndl,'XData',   p_foot1(1)                 , 'YData', p_foot1(2));
        set(Tibia1Hndl,'XData', [p_foot1(1) p_knee1(1)]   , 'YData', [p_foot1(2) p_knee1(2)]);
        set(Knee1Hndl,'XData',   p_knee1(1)               , 'YData', p_knee1(2));
        set(Femur1Hndl,'XData', [p_knee1(1) p_hip(1)]   , 'YData', [p_knee1(2) p_hip(2)]);
        set(Hip_Hndl,'XData',    p_hip(1)                 , 'YData', p_hip(2));
        set(Torso_Hndl,'XData', [p_hip(1) p_head(1)]    , 'YData', [p_hip(2) p_head(2)]);
        set(Femur2Hndl,'XData', [p_hip(1) p_knee2(1)]   , 'YData', [p_hip(2) p_knee2(2)]);
        set(Knee2Hndl,'XData',   p_knee2(1)               , 'YData', p_knee2(2));
        set(Tibia2Hndl,'XData', [p_knee2(1) p_foot2(1)] , 'YData', [p_knee2(2) p_foot2(2)]);
        set(foot2Hndl,'XData',   p_foot2(1)               , 'YData',  p_foot2(2));

        if(time(i)==T_impact(ImpactCounter))
            
            p_foot1=p_foot2;
            ImpactCounter=ImpactCounter+1;
        
            if(p_foot1(1)>XLimits(2))
                XLimits=XLimits+11;
                axis([get(gca, 'XLim')+11 -ViewWin ViewWin]);
            end
                
            
            if(mod(ImpactCounter,2)==0)

                uistack(Femur1Hndl, 'top')
                uistack(Knee1Hndl, 'top')
                uistack(Tibia1Hndl, 'top')
                uistack(foot1Hndl, 'top')
                set(Femur1Hndl,'color', 'g');
                set(Knee1Hndl,'MarkerFaceColor','g','MarkerEdgeColor','g');
                set(Tibia1Hndl,'color', 'g');
                set(foot1Hndl,'MarkerFaceColor','g','MarkerEdgeColor','g');

                set(Femur2Hndl,'color', 'b');
                set(Knee2Hndl,'MarkerFaceColor','b','MarkerEdgeColor','b');
                set(Tibia2Hndl,'color', 'b');
                set(foot2Hndl,'MarkerFaceColor','b','MarkerEdgeColor','b');

            else

                uistack(Femur2Hndl, 'top')
                uistack(Knee2Hndl, 'top')
                uistack(Tibia2Hndl, 'top')
                uistack(foot2Hndl, 'top')

                set(Femur1Hndl,'color', 'b');
                set(Knee1Hndl,'MarkerFaceColor','b','MarkerEdgeColor','b');
                set(Tibia1Hndl,'color', 'b');
                set(foot1Hndl,'MarkerFaceColor','b','MarkerEdgeColor','b');

                set(Femur2Hndl,'color', 'g');
                set(Knee2Hndl,'MarkerFaceColor','g','MarkerEdgeColor','g');
                set(Tibia2Hndl,'color', 'g');
                set(foot2Hndl,'MarkerFaceColor','g','MarkerEdgeColor','g');
            end
        end
        
        drawnow;
        pause((time(i+1)-time(i))*PaTi)
        
%         if(mod(i,10)==1)
%             frame=frame+1;
%             SavePdfFast(sprintf('%s/Walk_%03d', FramesFolder,frame))
%         end
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