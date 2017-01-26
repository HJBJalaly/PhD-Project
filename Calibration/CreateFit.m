function [fitresult, gof] = CreateFit(Q1, Q2, Value,FitType,DataType,fShow)

%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( Q1, Q2, Value );

% Set up fittype and options.
ft = fittype( FitType );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
opts.Robust = 'Bisquare';
opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
if(nargin==6)
    if(strcmp(fShow,'ShowOn'))
        figure( 'Name',['Type: ',FitType]);
        h = plot( fitresult, [xData, yData], zData );
        legend( h, 'untitled fit 1', [DataType, ' vs. Q1, Q2'], 'Location', 'NorthEast');
        % Label axes
        xlabel( 'Q1' );
        ylabel( 'Q2' );
        zlabel( DataType );
        grid on
        view( -118.5, 34.0 );
    end
end

