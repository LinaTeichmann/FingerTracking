function plotGeneralMeasures(gm,figuresDir)
%%


% Select variables
gm_names = fieldnames(gm);
y_labels = {'Movement time (sec.)', 'Endpoint Error', 'Endpoint Bias'};

% Organize colors
load('cdcol.mat')
C(1, 1, :) = cdcol.orange;
C(2, 1, :) = cdcol.cobaltblue;
C(3, 1, :) = cdcol.orange+0.07;
C(4, 1, :) = cdcol.lightblue;
C(5, 1, :) = cdcol.emeraldgreen;
C(6, 1, :) = cdcol.yellowgeen;


%% Plot
figureDim = [0 0 1 .5];
figure('units','normalized','outerposition',figureDim)

count = 1;
for i = [1 3 5]
    subplot(1,3,count)
    superbar(gm.(gm_names{i}), 'BarFaceColor', C, 'E', gm.(gm_names{i+1}), 'ErrorbarColor', C,'ErrorbarStyle', '|');
    set(gca, 'XLim', [.5 6.5]);
    
    if i == 1
%         set(gca, 'YLim', [max(gm.(gm_names{i}))*.8 max(gm.(gm_names{i})+gm.(gm_names{i+1}))], 'YTick', 0:0.1:max(gm.(gm_names{i})+gm.(gm_names{i+1})));
        set(gca, 'YLim', [max(gm.(gm_names{i}))*.8 1.005*max(gm.(gm_names{i})+gm.(gm_names{i+1}))]);
    elseif i == 5
        set(gca, 'YLim', [-0.7 0.7]);
    else
    end
    set(gca, 'XTickLabel', {'(A+B)+C','(AxB)+C','A+(B+C)','A+(BxC)','(A+B)-C','(AxB)-C'})
    set(gca, 'XTickLabelRotation', 25);
    set(gca, 'FontSize', 16);
    set(gcf, 'Color', 'white');
    ylabel(y_labels{count})
    box on;
    count = count+1;
end
savePNG(gcf,200, [figuresDir 'general_measures.png'])
close all

end

