function AnimBot3DOFOverplot(time,q,T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,PaTi)

FigName= 'One Step OverPlot';
% FramesFolder = './ImageExmaple';


% Initialize Figure
% -----------------

% initialize animation figure
AnimBot_Init(FigName);

% store figure handle for repeated access
FigHndl = findobj('Type', 'figure',  'Name', FigName);

% set figure axis range
ViewWin = 8*Lc_fem;  %[m]
axis([-.5 2 -.2 ViewWin]);


% Initialize Plot Handle, i.e. create plot dummy
% ----------------------------------------------


% Annotation: the Simulink inputs u(i) are not present at
% flag = 0. Hence, the plot dummy must be initiated with
% all values set to an arbitrary value (zero), i.e. a
% dummy is created.
d = [0; 0];

% MOHndl = plot(d, d, 'Color', [0.4 0.4 1], 'EraseMode','none', 'LineWidth', 12);
% GRHndl = plot(d, d, 'Color', 0.3*[1 1 1], 'EraseMode','none', 'LineWidth', 2);
% % foot 1
% foot1Hndl = plot(d, d, 'b', ...
%                     'MarkerFaceColor','b','MarkerEdgeColor','b',...
%                     'MarkerSize',10,...
%                     'Marker','o',...
%                     'LineStyle','none');
% % tibia 1
% Tibia1Hndl = plot(d, d, 'b',  'LineWidth',2);
% % knee 1
% Knee1Hndl = plot(d, d, 'b',  ...
%                     'MarkerFaceColor','b','MarkerEdgeColor','b',...
%                     'MarkerSize',10,...
%                     'Marker','o',...
%                     'LineStyle','none');
% % femur 1 
% Femur1Hndl = plot(d, d, 'b', 'LineWidth',2);
% % hip 
% Hip_Hndl = plot(d, d, 'r',  ...
%                     'MarkerFaceColor','r','MarkerEdgeColor','r',...
%                     'MarkerSize',10,...
%                     'Marker','o',...
%                     'LineStyle','none');
% % Torso
% Torso_Hndl = plot(d, d, 'y',  'LineWidth',2);
% % femur 2
% Femur2Hndl = plot(d, d, 'g',  'LineWidth',2);
% 
% % knee 2
% Knee2Hndl = plot(d, d, 'g', ...
%                     'MarkerFaceColor','g','MarkerEdgeColor','g',...
%                     'MarkerSize',10,...
%                     'Marker','o',...
%                     'LineStyle','none');
% 
% % tibia 2
% Tibia2Hndl = plot(d, d, 'g', 'LineWidth',2);
% % foot 2
% foot2Hndl = plot(d, d, 'g',  ...
%                     'MarkerFaceColor','g','MarkerEdgeColor','g',...
%                     'MarkerSize',10,...
%                     'Marker','o',...
%                     'LineStyle','none');
%                 
%                 
% TrackHandl= plot(d, d, 'Color',[1 1 1], 'EraseMode','none', ...
%     'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1],...
%     'MarkerSize',1,...
%     'Marker','none',...
%     'LineStyle','-.');
% 
% % link handles to current figure and to each other
% set(gca,   'UserData', foot1Hndl);
% set(foot1Hndl,'UserData', Tibia1Hndl);
% set(Tibia1Hndl,'UserData', Knee1Hndl);
% set(Knee1Hndl,'UserData', Femur1Hndl);
% set(Femur1Hndl,'UserData', Hip_Hndl);
% set(Hip_Hndl,'UserData', Torso_Hndl);
% set(Torso_Hndl,'UserData', Femur2Hndl);
% set(Femur2Hndl,'UserData', Knee2Hndl);
% set(Knee2Hndl,'UserData', Tibia2Hndl);
% set(Tibia2Hndl,'UserData', foot2Hndl);
% set(foot2Hndl,'UserData', GRHndl);
% set(GRHndl,'UserData', MOHndl);
% set(MOHndl,'UserData', TrackHandl);


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
    

    hold on
    Tracker=[];
    ImpactCounter=1;
    p_foot1    =[0 0]';
    LinsT1='-';
    LinsT2='-.';
    
    for j=4:-1:1
        IndXX1=find(time==T_impact(end-j),1,'last');
        IndXX2=find(time==T_impact(end-j+1),1,'last');

        for i=round(linspace(IndXX1,IndXX2-1,5))

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

            plot(   [XLimits(1) XLimits(2)]    , [0 0],'color',[0 0 0]);
    %         set(foot1Hndl,'XData',   p_foot1(1)                 , 'YData', p_foot1(2));
            plot([p_foot1(1) p_knee1(1)] , [p_foot1(2) p_knee1(2)],LinsT1,'linewidth',2);
    %         set(Knee1Hndl,'XData',   p_knee1(1)               , 'YData', p_knee1(2));
            plot([p_knee1(1) p_hip(1)]   , [p_knee1(2) p_hip(2)],LinsT1,'linewidth',2);
    %         set(Hip_Hndl,'XData',    p_hip(1)                 , 'YData', p_hip(2));
            plot([p_hip(1) p_head(1)] , [p_hip(2) p_head(2)],'linewidth',2);
            plot( [p_hip(1) p_knee2(1)]  , [p_hip(2) p_knee2(2)],LinsT2,'linewidth',2);
    %         set(Knee2Hndl,'XData',   p_knee2(1)               , 'YData', p_knee2(2));
            plot( [p_knee2(1) p_foot2(1)] , [p_knee2(2) p_foot2(2)],LinsT2,'linewidth',2);

            if(time(i)==T_impact(end-j+1))

                p_foot1L=p_foot1;
                p_foot1=p_foot2;
                LinsTemp=LinsT1;
                LinsT1=LinsT2;
                LinsT2=LinsTemp;
                
                if(mod(ImpactCounter,2)==0)

                else

                end
            end

            drawnow;
        end
        
        
        if(j==4 || j==2)
            i=IndXX1:IndXX2-1;
            q_fem1  = q(i,1)+q(i,5);
            q_fem2  = q(i,2)+q(i,5);
            q_tib1  = q(i,1)+q(i,3)+q(i,5);
            q_tib2  = q(i,2)+q(i,4)+q(i,5);

            p_knee1L    =repmat(p_foot1L,1,length(i)) + (L_tib          )*[ sin(pi-q_tib1), cos(pi-q_tib1)]';
            p_hipL      =p_knee1L+ (L_fem           )*[ sin(pi-q_fem1), cos(pi-q_fem1)]';
            p_knee2L    =p_hipL  + (L_fem           )*[ sin(q_fem2-pi),-cos(q_fem2-pi)]';
            p_foot2L    =p_knee2L+ (L_tib           )*[ sin(q_tib2-pi),-cos(pi-q_tib2)]';

            plot( p_foot2L(1,:) , p_foot2L(2,:),'color',[.8 .8 .8],'linewidth',2);
        end
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
        'Color',                        [1 1 1], ...
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
    'Color',     [1 1 1], ...
    'XColor',    [0 0 0], ...
    'YColor',    [0 0 0]);

axis on;
axis image;

hold on;
end