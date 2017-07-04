function xbyt_explore(expData, all, figuresDir)
% ar.calc1.generalMeasuresPlot('OM_exp2/sub','OM_exp2/add', 2, 0);
% sepAdd = 0 means no separation between addLS vs addSL
%

%%

%  ar.triplecode.xbyt_explore(AtB_pC)
%  subinitials = {'ac'}
%  experiment = 2


%% Preparing the data
subjIDs = tt.inf.listInitials(expData.d);
if all == 1
    subjIDs = {'all'};
end

originalTargets = 1:19;
targetXCoords = tt.nl.numberToX(originalTargets, 20);

% Condition
CondName = expData.d.all.Trials(1).Custom.Type;


%% Plotting
trajcolor = flipud(parula(19));
lineWidth = 2;
lineWidthDashed = 1;
textYCoord = 2.5;

for e = 1:length(subjIDs)
    hFig = figure(e);
    set(hFig, 'Position', [0 0 900 600])
    subData = expData.d.(subjIDs{e});
    targets = unique(arrayfun(@(t)t.Custom.ResTarget, subData.Trials));
    for i=1:length(targets);
        relevantTrials = subData.Trials(arrayfun(@(t)t.Custom.ResTarget==targets(i), subData.Trials));
        traj{i} = tt.preprocess.averageTrajectoryByTime(relevantTrials, .01, 0, [TrajCols.X, TrajCols.Y], 'mean');
        
        TrajColX = 2;
        
        TrajColY = 3;
        TrajColtime = 1;
        
        hold on
        plot(traj{i}(:,TrajColX), traj{i}(:,TrajColtime), 'Color', trajcolor(targets(i),:), 'LineWidth', lineWidth);
        connX = [traj{i}(end,TrajColX), targetXCoords(targets(i))];
        connY = [traj{i}(end,TrajColtime) textYCoord];
        plot(connX, connY, 'Color', trajcolor(targets(i),:), 'LineStyle', ':', 'LineWidth', lineWidthDashed);
        % Store movement time
        %             movtime{e}(i) = traj{i}(end,TrajColtime);
        
    end
    
    for itextTarg= 1:length(originalTargets)
        text(targetXCoords(itextTarg), textYCoord, sprintf('%2d',originalTargets(itextTarg)), 'Color', trajcolor(itextTarg,:), 'HorizontalAlignment', 'center', 'FontSize', 24);
    end
    
    % Title, axes, gridlines
    set(gca, 'YTick', 0:0.2:textYCoord-0.2);
    nTicks = floor(length(1:9)/10);
    tickSize = length(1:9) * 2 / nTicks;
    xTicks = -length(1:19) : tickSize : length(1:19);
    set(gca, 'XLim', [-.75 .75]);
    set(gca, 'YLim', [0 textYCoord+0.1]);
    set(gca, 'XTick', xTicks);
    set(gca, 'XTickLabel', arrayfun(@(n){sprintf('%d', (n-1)*10)}, 1:(nTicks+1)));
    set(gca, 'TickLength', [0 0]);
    set(gca, 'FontSize', 30);
    set(gcf, 'Color', 'white');
    xlabh= xlabel('X coordinate');
    set(xlabh,'Position',get(xlabh,'Position') - [0 0.02 0])
    ylabel('Time (sec.)')
    grid on
    box on;
    title([CondName '  subject  ' subjIDs{e}])
    savePNG(gcf,200, [figuresDir CondName '_' subjIDs{e} '.png'])
    close all
end
end


