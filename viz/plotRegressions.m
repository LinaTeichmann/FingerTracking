function plotRegressions(rr,ds_name,reg_name, figuresDir, varargin)

cdcol = load('cdcol.mat');
cdcol = cdcol.cdcol;

%% Get arguments varargin
[y_lim, max_time] = parseArgs(varargin);

% Plot
if isempty(y_lim); y_lim = [-.35 1.1]; else end
if isempty(max_time); max_time = 2.2; else end

%% Colors


%% Prepare predictors
% Title
ds_title = ds_name{2};
ds_name = ds_name{1};
% Add .b in each predictor
for i = 1:length(rr.(ds_name).avg.(reg_name).predictorNames);
    predictors{i} = [rr.(ds_name).avg.(reg_name).predictorNames{i} '.b'];
end
predictors_plot = predictors(2:end-2);

%% Define colors and other parameters
% Colors list
color_list = {cdcol.ultramarine,cdcol.grassgreen,cdcol.orange,cdcol.scarlet, cdcol.lightblue};
colors_plot = color_list(1:length(predictors_plot));
colors_plot(end+1) = {[0 0 0]};
colors_plot(end+1) = {[.5 .5 .5]};

%% Plot
cmp = tt.reg.compareParams(rr.(ds_name), reg_name, predictors_plot, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr.(ds_name), reg_name, predictors(1), 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', y_lim, 'print','basic','MaxTime', max_time, 'Color', colors_plot, 'marker', 'o','markersize', 12, 'nonsmarkers', 'noclf');
updatePlotParm(rr.(ds_name).avg.(reg_name).PredictorDesc(2:end-2), ds_title, max_time)

%% Save
savePNG(gcf,200, [figuresDir, 'reg_', ds_name, '_' reg_name,'.png'])
%close all



%% Parse args
    function [y_lim, max_time] = parseArgs(args)
        % Prelocate args
        y_lim = [];
        max_time = [];
        args = stripArgs(args);
        while ~isempty(args)
            switch(lower(args{1}))
                case 'y_lim'
                    y_lim = args{2};
                    args = args(2:end);
                case 'max_time'
                    max_time = args{2};
                    args = args(2:end);
                otherwise
                    error('Unsupported argument "%s"!', args{1});
            end
            args = stripArgs(args(2:end));
        end
    end


end

