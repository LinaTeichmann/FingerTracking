%% Finger-tracking MEG  

%% Initialize Directiories
addpath(genpath('/Users/pinheirochagas/Pedro/NeuroSpin/Experiments/iPad_Calc/trajtracker/matlab/'))
addpath(genpath('/Users/pinheirochagas/Pedro/NeuroSpin/Experiments/iPad_Calc/scripts/Calc_iPad/'))
addpath(genpath('/Users/pinheirochagas/Pedro/NeuroSpin/Experiments/Calculia/scripts/Calc_ECoG/PedroMatlabCustom/'))
%% Figures path
figuresDir = [TrajTrackerDataPath '/arit_express/results/figures/'];


%% 1. Preprocessing
tt.preprocessSet('arit_express/data/arit_express/two_operations/all', 'ProcessDSFunc', @postprocess_set)

% Load data sets
addsub = tt.loadDataset('arit_express/data/arit_express/two_operations/all');

%% Dynamic filtering
add = tt.util.filterDataset(addsub, @(t) strcmp(t.Custom.Type,'A+B') == 1, 'condName', 'A+B'); 
sub = tt.util.filterDataset(addsub, @(t) strcmp(t.Custom.Type,'A-B') == 1, 'condName', 'A-B'); 

%% General plot
for i = 1:length(ds_names)
    % Individual data
    xbyt_explore(ds.(ds_names{i}), 0, [figuresDir 'trajectories/individual/']);
    % Avarage across subjects
    xbyt_explore(ds.(ds_names{i}), 1, [figuresDir 'trajectories/group/']);
end


%% Geneneral Performance Measures
% generalPerfc1 = tt.inf.printBasicStats({c1addSL,c1addLS,c1sub,c11d});
generalPerf = tt.inf.printBasicStats({add,sub,addsub});

% Get mean and std
for i = 1:length(generalPerf)
    gm.MTMean(i) = mean(generalPerf{i}.MovementTime);
    gm.MTSD(i) = std(generalPerf{i}.MovementTime)/sqrt(nsubjects);
    gm.EPEMean(i) = mean(generalPerf{i}.EndpointErr);
    gm.EPESD(i) = std(generalPerf{i}.EndpointErr)/sqrt(nsubjects);
    gm.EPBMean(i) = mean(generalPerf{i}.EndpointBias);
    gm.EPBSD(i) = std(generalPerf{i}.EndpointBias)/sqrt(nsubjects);
end

%% Plotting general measures
plotGeneralMeasures(gm,figuresDir)



%% Regressions
% Define regression parameters
regSampRate = 0.03;
regMaxTime = 2.5;
% Main regressions
for i = 1:length(ds_names);
    rr1.(ds_names{i}) = runRegressions(ds.(ds_names{i}), 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate}, 'saveas', ['rr1_' ds_names{i}, '.mat']);
end

%% Load regressions
for i=1:length(ds_names)
    reg_tmp = load([TrajTrackerDataPath '/arit_express/results/regressions/' 'rr1_' ds_names{i,1} '.mat']);
    rr1.(ds_names{i,1}) = reg_tmp.rr;
end

%% Plot regressions
for i = 1:length(ds_names)
    get_regs = fieldnames(rr1.(ds_names{i}).avg);
    get_regs = get_regs(5:end);
    for ii = 1:length(get_regs)
        plotRegressions(rr1,ds_names(i,:),get_regs{ii}, figuresDir)
    end
end







%% Update velocity (include implied endpoint velocity)
twoOps = updateVelAcc(twoOps, 'iep', 0.02);
add_0_20 = updateVelAcc(add_0_20, 'iep', 0.02);
mult_0_20 = updateVelAcc(mult_0_20, 'iep', 0.02);

%% Add acceleration bursts
% tt.vel.findAccelerationBursts(twoOps.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'MinTime', 0.1, 'Axis', 'Y', 'AttrPrefix', 'Y')
% tt.vel.findAccelerationBursts(twoOps.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'MinTime', 0.1, 'Axis', 'X', 'AttrPrefix', 'X')
tt.vel.findAccelerationBursts(twoOps.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'Y', 'AttrPrefix', 'Y')
tt.vel.findAccelerationBursts(twoOps.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'X', 'AttrPrefix', 'X')
tt.vel.findAccelerationBursts(twoOps.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'iep', 'AttrPrefix', 'X')

tt.vel.findAccelerationBursts(add_0_20.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'Y', 'AttrPrefix', 'Y')
tt.vel.findAccelerationBursts(add_0_20.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'X', 'AttrPrefix', 'X')
tt.vel.findAccelerationBursts(add_0_20.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'iep', 'AttrPrefix', 'X')

tt.vel.findAccelerationBursts(mult_0_20.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'Y', 'AttrPrefix', 'Y')
tt.vel.findAccelerationBursts(mult_0_20.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'X', 'AttrPrefix', 'X')
tt.vel.findAccelerationBursts(mult_0_20.d, 'MinDur', 0.07, 'MinAcc', 0.2, 'Smooth', 0.1, 'Axis', 'iep', 'AttrPrefix', 'X')

% Update acceleration burst info
filterAccBurst(twoOps.d)
filterAccBurst(add_0_20.d)
filterAccBurst(mult_0_20.d)



%% Dynamic filtering
ds.ApB_pC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'(A+B)+C') == 1, 'condName', '(A+B)+C'); 
ds.AtB_pC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'(AxB)+C') == 1, 'condName', '(AxB)+C'); 
ds.Ap_BpC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'A+(B+C)') == 1, 'condName', 'A+(B+C)'); 
ds.Ap_BtC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'A+(BxC)') == 1, 'condName', 'A+(BxC)'); 
ds.ApB_mC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'(A+B)-C') == 1, 'condName', '(A+B)-C'); 
ds.AtB_mC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'(AxB)-C') == 1, 'condName', '(AxB)-C'); 
ds_names = fieldnames(ds);
ds_names(1,2) = {'(A+B)+C'}; ds_names(2,2) = {'(AxB)+C'}; ds_names(3,2) = {'A+(B+C)'}; ds_names(4,2) = {'A+(BxC)'}; ds_names(5,2) = {'(A+B)-C'}; ds_names(6,2) = {'(AxB)-C'};
nsubjects = length(tt.inf.listInitials(twoOps.raw));

ApB_pC_ApB_mC = tt.util.filterDataset(twoOps, @(t) strcmp(t.Custom.Type,'(A+B)+C') == 1 | strcmp(t.Custom.Type,'(A+B)-C') == 1, 'condName', '(A+B)+/-C'); 


%% General plot
for i = 1:length(ds_names)
    % Individual data
    xbyt_explore(ds.(ds_names{i}), 0, [figuresDir 'trajectories/individual/']);
    % Avarage across subjects
    xbyt_explore(ds.(ds_names{i}), 1, [figuresDir 'trajectories/group/']);
end

%% Plot acceleration
for i = 1:length(ds_names)
    % Individual data
    plotFilledTraj(ds.(ds_names{i}), 0, [figuresDir 'trajectories/individual/']);
end

%% Geneneral Performance Measures
% generalPerfc1 = tt.inf.printBasicStats({c1addSL,c1addLS,c1sub,c11d});
generalPerf = tt.inf.printBasicStats({ds.ApB_pC,ds.AtB_pC,ds.Ap_BpC,ds.Ap_BtC,ds.ApB_mC,ds.AtB_mC});

% Get mean and std
for i = 1:length(generalPerf)
    gm.MTMean(i) = mean(generalPerf{i}.MovementTime);
    gm.MTSD(i) = std(generalPerf{i}.MovementTime)/sqrt(nsubjects);
    gm.EPEMean(i) = mean(generalPerf{i}.EndpointErr);
    gm.EPESD(i) = std(generalPerf{i}.EndpointErr)/sqrt(nsubjects);
    gm.EPBMean(i) = mean(generalPerf{i}.EndpointBias);
    gm.EPBSD(i) = std(generalPerf{i}.EndpointBias)/sqrt(nsubjects);
end

%% Plotting general measures
plotGeneralMeasures(gm,figuresDir)


%% Plotting
for i = 1:length(ds_names)
    get_regs = fieldnames(rr1.(ds_names{i}).avg);
    get_regs = get_regs(5:end);
    for ii = 1:length(get_regs)
        plotRegressions(rr1,ds_names(i,:),get_regs{ii}, figuresDir)
    end
end

%% Regressions
% Define regression parameters
regSampRate = 0.03;
regMaxTime = 2.5;
% Main regressions
for i = 1:length(ds_names);
    rr1.(ds_names{i}) = runRegressions(ds.(ds_names{i}), 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate}, 'saveas', ['rr1_' ds_names{i}, '.mat']);
end
% 
% for i = 1:length(ds_names);
%     rr3.(ds_names{i}) = runRegressions(ds.(ds_names{i}), 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate}, 'saveas', [ds_names{i}, '.mat']);
% end
% 
% for i = 1:length(ds_names);
%     rr_alig.(ds_names{i}) = runRegressions(ds.(ds_names{i}), 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate, 'Row1OffsetFunc', @(trial) trial.Custom.YPosAccelStartRows(1), 'TrialFilter', @(trial) ~isempty(trial.Custom.YPosAccelStartRows)}, 'saveas', [ds_names{i}, '_alig_acc_burst.mat']);
% end

% ApB_pC_ApB_mC
rr_comb.('ApB_pC_ApB_mC') = runRegressions(ApB_pC_ApB_mC, 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate}, 'saveas', ['rr1_' ds_names{i}, '.mat']);
rr_comb2.('ApB_pC_ApB_mC') = runRegressions(ApB_pC_ApB_mC, 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate}, 'saveas', ['rr1_' ds_names{i}, '.mat']);

%% Load regressions
for i=1:length(ds_names)
    reg_tmp = load([TrajTrackerDataPath '/arit_express/results/regressions/' 'rr1_' ds_names{i,1} '.mat']);
    rr1.(ds_names{i,1}) = reg_tmp.rr;
end


%% Plot regressions
for i = 1:length(ds_names)
    get_regs = fieldnames(rr1.(ds_names{i}).avg);
    get_regs = get_regs(5:end);
    for ii = 1:length(get_regs)
        plotRegressions(rr1,ds_names(i,:),get_regs{ii}, figuresDir)
    end
end



%% Plot some individual regressions
plotRegressions(rr1,ds_names(6,:),get_regs{4},figuresDir, 'y_lim', [-1.5 5.5])

plotRegressions(rr_comb,{'ApB_pC_ApB_mC', '(A+B)+/-C'},'Op12sign3abs3',figuresDir)

plotRegressions(rr_comb2,{'ApB_pC_ApB_mC', '(A+B)+/-C'},'Op12sign3abs3operator2',figuresDir)


%% Get some timing of the effects
beta_val = 0.5;

rr1.AtB_pC.avg.OutPar.times(find(rr1.AtB_pC.avg.OutPar.getParamValue('ResParent', 'b') >= beta_val, 1))
rr1.AtB_pC.avg.OutPar.times(find(rr1.AtB_pC.avg.OutPar.getParamValue('OpOutParent', 'b') >= beta_val, 1))
rr1.AtB_pC.avg.OutPar.times(find(rr1.AtB_pC.avg.OutPar.getParamValue('OpOutParent', 'b') >= beta_val, 1)) - rr1.AtB_pC.avg.OutPar.times(find(rr1.AtB_pC.avg.OutPar.getParamValue('ResParent', 'b') >= beta_val, 1)) 


rr1.ApB_pC.avg.OutPar.times(find(rr1.ApB_pC.avg.OutPar.getParamValue('ResParent', 'b') >= beta_val, 1))
rr1.ApB_pC.avg.OutPar.times(find(rr1.ApB_pC.avg.OutPar.getParamValue('OpOutParent', 'b') >= beta_val, 1))
rr1.ApB_pC.avg.OutPar.times(find(rr1.ApB_pC.avg.OutPar.getParamValue('OpOutParent', 'b') >= beta_val, 1)) - rr1.ApB_pC.avg.OutPar.times(find(rr1.ApB_pC.avg.OutPar.getParamValue('ResParent', 'b') >= beta_val, 1))


rr1.ApB_mC.avg.OutPar.times(find(rr1.ApB_mC.avg.OutPar.getParamValue('ResParent', 'b') >= beta_val, 1))
rr1.ApB_mC.avg.OutPar.times(find(rr1.ApB_mC.avg.OutPar.getParamValue('OpOutParent', 'b') >= beta_val, 1))
rr1.ApB_mC.avg.OutPar.times(find(rr1.ApB_mC.avg.OutPar.getParamValue('OpOutParent', 'b') >= beta_val, 1)) - rr1.ApB_mC.avg.OutPar.times(find(rr1.ApB_mC.avg.OutPar.getParamValue('ResParent', 'b') >= beta_val, 1))




round(mean(acc.AtB_pC.yacc_nbst2_pos_time))
round(mean(acc.AtB_pC.yacc_nbst3_pos_time))
round(mean(acc.AtB_pC.yacc_nbst3_pos_time)) - round(mean(acc.AtB_pC.yacc_nbst2_pos_time))

round(mean(acc.ApB_pC.yacc_nbst2_pos_time))
round(mean(acc.ApB_pC.yacc_nbst3_pos_time))
round(mean(acc.ApB_pC.yacc_nbst3_pos_time)) - round(mean(acc.ApB_pC.yacc_nbst2_pos_time))

round(mean(acc.ApB_mC.yacc_nbst2_pos_time))
round(mean(acc.ApB_mC.yacc_nbst3_pos_time))
round(mean(acc.ApB_mC.yacc_nbst3_pos_time)) - round(mean(acc.ApB_pC.yacc_nbst2_pos_time))


%% Compare (A+B)+C and (AxB)+C

figure(1)
cmp = tt.reg.compareParams({rr1.ApB_pC, rr1.AtB_pC}, 'OutPar', {'OpOutParent.b','OpOutParent.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2.2, 'Color', {cdcol.yellow,cdcol.orange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 12, 'nonsmarkers');
updatePlotParm({'(A+B)+C', '(AxB)+C'}, 'out()')
savePNG(gcf,200, [figuresDir, '(A+B)+C_(AxB)+C_Out', '.png'])

figure(2)
cmp = tt.reg.compareParams({rr1.ApB_pC, rr1.AtB_pC}, 'OutPar', {'ResParent.b','ResParent.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2.2, 'Color', {cdcol.yellow,cdcol.orange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 12, 'nonsmarkers');
updatePlotParm({'(A+B)+C', '(AxB)+C'}, 'result()')
savePNG(gcf,200, [figuresDir, '(A+B)+C_(AxB)+C_Par', '.png'])

%% Compare A+(B+C) A+(BxC)

figure(1)
cmp = tt.reg.compareParams({rr1.Ap_BpC, rr1.Ap_BtC}, 'OutPar', {'OpOutParent.b','OpOutParent.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2.2, 'Color', {cdcol.yellow,cdcol.orange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 12, 'nonsmarkers');
updatePlotParm({'A+(B+C)', 'A+(BxC)'}, 'out()')
savePNG(gcf,200, [figuresDir, 'A+(B+C)_A+(BxC)_Out', '.png'])

figure(2)
cmp = tt.reg.compareParams({rr1.Ap_BpC, rr1.Ap_BtC}, 'OutPar', {'ResParent.b','ResParent.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2.2, 'Color', {cdcol.yellow,cdcol.orange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 12, 'nonsmarkers');
updatePlotParm({'A+(B+C)', 'A+(BxC)'}, 'result()')
savePNG(gcf,200, [figuresDir, 'A+(B+C)_A+(BxC)_Par', '.png'])

%% plot velocity profiles
hold on
tt.vis.plotTrajValue(ds.ApB_pC.d, TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(ds.AtB_pC.d, 'TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'Colors', {'k'}, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))

hold on
tt.vis.plotTrajValue(ds.Ap_BpC.d, 'TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(ds.Ap_BtC.d, 'TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'Colors', {'k'}, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))

hold on
tt.vis.plotTrajValue(ds.Ap_BpC.d, 'TrajCol', TrajCols.AngularVelocity, 'GrpFunc', @(trial)1, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(ds.Ap_BtC.d, 'TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'Colors', {'k'}, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))


hold on
tt.vis.plotTrajValue(ds.ApB_pC.d, 'GetValueFunc', @getAngularAcc, 'GrpFunc', @(trial)1, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(ds.AtB_pC.d, 'GetValueFunc', @getAngularAcc, 'GrpFunc', @(trial)1, 'Colors', {'k'}, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))


hold on
tt.vis.plotTrajValue(ds.ApB_pC.d, 'TrajCol', TrajCols.YVelocity, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(ds.AtB_pC.d, 'TrajCol', TrajCols.YVelocity, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                         Add Mult                        %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analyze add_mult to compare the movement time
targets_add_0_20 = unique(arrayfun(@(t){t.Custom.PresentedTarget}, add_0_20.d.all.Trials));
results_add_0_20 = unique(arrayfun(@(t)t.Custom.ResTarget, add_0_20.d.all.Trials));

targets_mult_0_20 = unique(arrayfun(@(t){t.Custom.PresentedTarget}, mult_0_20.d.all.Trials));
results_mult_0_20 = unique(arrayfun(@(t)t.Custom.ResTarget, mult_0_20.d.all.Trials));
results_mult_0_20 = results_mult_0_20(1:end-1);

count = 1;
for i = results_mult_0_20
    MT_add_0_20(:,count) = tt.inf.getTrialValue(add_0_20.d, 'Getter',  @(trial) trial.MovementTime, 'PerSubj', 'TrialFilter', @(trial) trial.Custom.ResTarget == i);
    MT_mult_0_20(:,count) = tt.inf.getTrialValue(mult_0_20.d, 'Getter', @(trial) trial.MovementTime, 'PerSubj', 'TrialFilter', @(trial) trial.Custom.ResTarget == i);
    count = count + 1;
end

MT_add_0_20_mean = mean(MT_add_0_20,1);
MT_add_0_20_std = std(MT_add_0_20,1)/sqrt(length(MT_add_0_20_mean));

MT_mult_0_20_mean = mean(MT_mult_0_20,1);
MT_mult_0_20_std = std(MT_mult_0_20,1)/sqrt(length(MT_mult_0_20_mean));

% Plot
figureDim = [0 0 .7 .55];
figure('units','normalized','outerposition',figureDim)
C = [cdcol.ultramarine; cdcol.grassgreen];
superbar(results_mult_0_20, [MT_add_0_20_mean;MT_mult_0_20_mean]', 'BarFaceColor', permute(C, [3 1 2]), 'E', [MT_add_0_20_std;MT_mult_0_20_std]', 'ErrorbarColor', permute(C, [3 1 2]),'ErrorbarStyle', '|');
set(gca, 'FontSize', 20)
ylim([1, 1.5])
ylabel('Movement Time (sec.)')
xlabel('Results')
% legend({'Additions', 'Multiplications'}, 'Location', 'SouthEast')
savePNG(gcf,200, [figuresDir, 'add_mult', '.png'])

%% Regressions
% Define regression parameters
am.add_0_20 = add_0_20;
am.mult_0_20 = mult_0_20;
am_names = fieldnames(am);
am_names{1,2} = 'A+B';
am_names{2,2} = 'AxB';

regSampRate = 0.03;
regMaxTime = 2;
% Main regressions
for i = 1:length(am_names);
    rr_addmult.(am_names{i}) = runRegressions(am.(am_names{i}), 'RegArgs', {'MaxTime', regMaxTime, 'dt', regSampRate}, 'saveas', ['rr_addmult_' am_names{i}, '.mat']);
end

%% Load regressions
for i=1:length(am_names)
    reg_tmp = load([TrajTrackerDataPath '/arit_express/results/regressions/' 'rr_addmult_' am_names{i,1} '.mat']);
    rr_addmult.(am_names{i,1}) = reg_tmp.rr;
end



%% Plot regressions
for i = 1:length(am_names)
    get_regs = fieldnames(rr_addmult.(am_names{i}).avg);
    get_regs = get_regs(5:end);
    for ii = 1:length(get_regs)
        plotRegressions(rr_addmult,am_names(i,:),get_regs{ii}, figuresDir)
    end
end


plotRegressions(rr_addmult,am_names(1,:),get_regs{1},figuresDir, 'y_lim', [-1.5 5.5])

plotRegressions(rr_addmult,am_names(2,:),get_regs{3},figuresDir, 'max_time', 1.5, 'y_lim', [-0.35 1.1])
close all

ylim([-1.2 5])

%% Compare result regression add and mult
figure(1)
cmp = tt.reg.compareParams({rr_addmult.add_0_20, rr_addmult.mult_0_20}, 'Target', {'ResTarget.b','ResTarget.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 1.5, 'Color', {cdcol.yellow,cdcol.orange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 12, 'nonsmarkers');
updatePlotParm({'A+B', 'AxB'}, 'Result', 1.5)
savePNG(gcf,200, [figuresDir, 'A+B_AxB_result', '.png'])

%% Compare result regression add and mult
figure(1)
cmp = tt.reg.compareParams({rr_addmult.add_0_20, rr_addmult.mult_0_20}, 'Target', {'ResTarget.beta','ResTarget.beta'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 1.5, 'Color', {cdcol.yellow,cdcol.orange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 12, 'nonsmarkers');
updatePlotParm({'A+B', 'AxB'}, 'Result', 1.5)
savePNG(gcf,200, [figuresDir, 'A+B_AxB_result', '.png'])

%% Y velocity as a function of time. 
hold on
tt.vis.plotTrajValue(add_0_20.d, 'TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(mult_0_20.d, 'TrajCol', TrajCols.YVelocity, 'GrpFunc', @(trial)1, 'Colors', {'k'}, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))


hold on
tt.vis.plotTrajValue(add_0_20.d, 'TrajCol', TrajCols.YAcceleration, 'GrpFunc', @(trial)1, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.vis.plotTrajValue(mult_0_20.d, 'TrajCol', TrajCols.YAcceleration, 'GrpFunc', @(trial)1, 'Colors', {'k'}, 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))



%% Get average movement time
tt.inf.getTrialValue(add_0_20.d, 'Attr', 'MovementTime', 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))
tt.inf.getTrialValue(mult_0_20.d, 'Attr', 'MovementTime', 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20))


add_mt = tt.inf.getTrialValue(add_0_20.d, 'Getter',  @(trial) trial.MovementTime, 'PerSubj', 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20));
mult_mt = tt.inf.getTrialValue(mult_0_20.d, 'Getter',  @(trial) trial.MovementTime, 'PerSubj', 'TrialFilter', @(trial)ismember(trial.Target, results_mult_0_20));


'TrialFilter', @(trial)ismember(trial.Target, [2 4 6])

    MT_mult_0_20(:,count) = tt.inf.getTrialValue(mult_0_20.d, 'Getter', @(trial) trial.MovementTime, 'PerSubj', 'TrialFilter', @(trial) trial.Custom.ResTarget == i);
    count = count + 1;
end




%%
figure(5)
cmp = tt.reg.compareParams(rr_Ap_BtC, 'Op123Parent', {'Operand1.b', 'Operand2.b', 'Operand3.b', 'ResParent.b' }, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_Ap_BtC, 'Op123Parent', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customGreen, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'A', 'B', 'C', 'AxB'}, 'A+(BxC)')


figure(4)
cmp = tt.reg.compareParams({rr_Ap_BtC, rr_AtB_pC}, 'Op123Parent', {'Operand1.b', 'Operand3.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers');
updatePlotParm({'A+(BxC)', '(AxB)+C'}, 'Outside Parentheses')

figure(5)
cmp = tt.reg.compareParams({rr_Ap_BtC, rr_AtB_pC}, 'Op123Parent', {'ResParent.b','ResParent.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers');
updatePlotParm({'A+(BxC)', '(AxB)+C'}, 'Parentheses')


cmp = tt.reg.compareParams({rr2.ApB_pC, rr2.AtB_pC}, 'OutPar', {'OpOutParent.b','OpOutParent.b'}, 'paramsMethod', 'Tzero');
cmp = tt.reg.compareParams({rr2.ApB_pC, rr2.AtB_pC}, 'OutPar', {'ResParent.b','ResParent.b'}, 'paramsMethod', 'Tzero');

tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {cdcol.turquoiseblue,cdcol.grassgreen, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers');





save2pdf([figuresDir 'rr_AtB_pC_Op1Op2Op3ParentPrevTldrfix.pdf'], gcf, 600)       




%%




Op123Parent
OutMinMaxParent


%%
figure(1)

cmp = tt.reg.compareParams(rr_AtB_pC, 'OutMinMaxParent', {'OpMinParent.b', 'OpMaxParent.b', 'OpOutParent.b', 'ResParent.b' }, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_AtB_pC, 'OutMinMaxParent', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customGreen, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'Min()', 'Max()', 'out()', '()result'}, '(AxB)+C')


%%
figure(2)
cmp = tt.reg.compareParams(rr_Ap_BtC, 'OutMinMaxParent', {'OpMinParent.b', 'OpMaxParent.b', 'OpOutParent.b', 'ResParent.b' }, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_Ap_BtC, 'OutMinMaxParent', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customGreen, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'Min()', 'Max()', 'out()', '()result'}, 'A+(BxC)')


figure(4)
cmp = tt.reg.compareParams({rr_Ap_BtC, rr_AtB_pC}, 'OutMinMaxParent', {'Operand1.b', 'Operand3.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers');
updatePlotParm({'A+(BxC)', '(AxB)+C'}, 'Outside Parentheses')

figure(5)
cmp = tt.reg.compareParams({rr_Ap_BtC, rr_AtB_pC}, 'OutMinMaxParent', {'ResParent.b','ResParent.b'}, 'paramsMethod', 'Tzero');
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers');
updatePlotParm({'A+(BxC)', '(AxB)+C'}, 'Parentheses')



%%
figure(1)
cmp = tt.reg.compareParams(rr_AtB_pC, 'OutSumParProdPar', {'OpOutParent.b','SumParent.b', 'ResParent.b'}, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_AtB_pC, 'OutSumParProdPar', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'out()', 'sum()', 'prod()'}, '(AxB)+C')




%%
figure(2)
cmp = tt.reg.compareParams(rr_Ap_BtC, 'OutSumParProdPar', {'OpOutParent.b','SumParent.b', 'ResParent.b'}, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_Ap_BtC, 'OutSumParProdPar', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'out()', 'sum()', 'prod()'}, 'A+(BxC)')



subplot(2,2,1)
figure(3)

cmp = tt.reg.compareParams(rr_AtB_mC, 'OutSumParProdPar', {'OpOutParent.b*-1','SumParent.b','ResParent.b' }, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_AtB_mC, 'OutSumParProdPar', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'out()', 'sum()', 'prod()'}, '(AxB)-C')






figure(1)

cmp = tt.reg.compareParams(rr_AtB_mC, 'Op123Parent', {'Operand1.b', 'Operand2.b', 'Operand3.b*-1', 'ResParent.b' }, 'paramsMethod', 'Tzero');
cmpConst = tt.reg.compareParams(rr_AtB_mC, 'Op123Parent', {'const.b'}, 'paramsMethod', 'Tzero');
plot(cmpConst.times,cmpConst.cmpParam.values/10,'Color',[.9 .9 .9], 'LineWidth',4, 'HandleVisibility', 'off')
tt.reg.plotParamComparison(cmp, 'legend','legend', 'StdErr', 'ShadeVar', .1, 'YLim', [-.35 1.1], 'print','basic','MaxTime', 2, 'Color', {customBlue,customOrange, customGreen, customPurple, [0 0 0], [.5 .5 .5]}, 'marker', 'o','markersize', 8, 'nonsmarkers', 'noclf');
updatePlotParm({'A', 'B', 'C', 'AxB'}, '(AxB)-C')




%% Averaged trajectories by result
ar.calc1.generalMeasuresPlot('OM_exp1/sub', 'OM_exp1/add',1,0);
figure(1)
% save2pdf([figuresDir 'xbyt_AddLS_AddSL_Sub_exp1.pdf'], gcf, 600)
save2pdf([figuresDir 'xbyt_Sub_ADD_exp1.pdf'], gcf, 600)

% Problem difficulty effect
ar.calc1.plotProblemDifficulty(1)
figure(2)
save2pdf([figuresDir 'ProblemDifficultyEffect_exp1.pdf'], gcf, 600)

% % Movement time, endpoint error and endpoint bias
% ar.calc1.plotMovTimeEPerrorEPbias(1) 
% save2pdf([figuresDir 'MovTimeEPErrorEPbias_exp1.pdf'], gcf, 600)











     

%% ApB_pC
cmp = treg.compareParams(rr_ApB_pC, 'ApB_pC_Op1Op2Op3PrevTldrfix', {'b_Operand1', 'b_Operand2', 'b_Operand3' }, 'compare1method', 'Tzero');
% cmp = treg.compareParams(rr_ApB_pC, 'ApB_pC_Op1Op2Op3PrevTldrfix', {'b_OpMinParent', 'b_OpMaxParent', 'b_OpOutParent'}, 'compare1method', 'Tzero');


% preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-1 1.25], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen});
preg.plotParamComparison(cmp, 'legend','legend','YLim', [-.25 1.2], 'print','basic', 'MaxTime', 2, 'Color', {customBlue,customOrange, customGreen});

ar.triplecode.updatePlotParm({'Operand 1', 'Operand 2', 'Operand 3'}, '(A+B)+C')
% ar.triplecode.updatePlotParm({'MaxOpParent', 'MinOpParent', 'OpOutParent'}, '(A+B)+C')
save2pdf([figuresDir 'rr_ApB_pC_Op1Op2Op3PrevTldrfix.pdf'], gcf, 600)  

% AtB_mC
cmp = treg.compareParams(rr_AtB_mC, 'AtB_mC_Op1Op2SignOp3ParentPrevTldrfix', {'b_Operand1', 'b_Operand2', 'b_SignOperand3', 'b_ResParent'}, 'compare1method', 'Tzero');
preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-1 1.6], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen, customCyan});
preg.plotParamComparison(cmp, 'legend','legend', 'YLim', [-.25 1.2], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen, customCyan});

ar.triplecode.updatePlotParm({'Operand 1', 'Operand 2', 'Operand 3', 'ResParent'}, '(AxB)-C')
save2pdf([figuresDir 'rr_AtB_mC_Op1Op2SignOp3ParentPrevTldrfix.pdf'], gcf, 600)  

% ApB_mC
cmp = treg.compareParams(rr_ApB_mC, 'ApB_mC_Op1Op2SignOp3PrevTldrfix', {'b_Operand1', 'b_Operand2', 'b_SignOperand3'}, 'compare1method', 'Tzero');
preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-.55 1.25], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen});
preg.plotParamComparison(cmp, 'legend','legend', 'YLim', [-.25 1.2], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen});


ar.triplecode.updatePlotParm({'Operand 1', 'Operand 2', 'Operand 3'}, '(A+B)-C')
save2pdf([figuresDir 'rr_AtB_pC_Op1Op2Op3ParentPrevTldrfix_dp.pdf'], gcf, 600)  

% Ap_BtC
cmp = treg.compareParams(rr_Ap_BtC, 'Ap_BtC_Op1Op2Op3ParentPrevTldrfix', {'b_Operand1', 'b_Operand2', 'b_Operand3', 'b_ResParent'}, 'compare1method', 'Tzero');
% cmp = treg.compareParams(rr_Ap_BtC, 'Ap_BtC_Op1Op2Op3ParentPrevTldrfix', {'b_OpMinParent', 'b_OpMaxParent', 'b_OpOutParent'}, 'compare1method', 'Tzero');

% preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-0.55 1.25], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen, customCyan});
preg.plotParamComparison(cmp, 'legend','legend', 'YLim', [-0.25 1.2], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen, customCyan});

ar.triplecode.updatePlotParm({'Operand 1', 'Operand 2', 'Operand 3', 'ResParent'}, 'A+(BxC)')
save2pdf([figuresDir 'rr_Ap_BtC_Op1Op2Op3ParentPrevTldrfix.pdf'], gcf, 600) 

% Ap_BpC
cmp = treg.compareParams(rr_Ap_BpC, 'Ap_BpC_Op1Op2Op3PrevTldrfix', {'b_Operand1', 'b_Operand2', 'b_Operand3'}, 'compare1method', 'Tzero');
% cmp = treg.compareParams(rr_Ap_BpC, 'Ap_BpC_Op1Op2Op3PrevTldrfix', {'b_OpMinParent', 'b_OpMaxParent', 'b_OpOutParent'}, 'compare1method', 'Tzero');

preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-0.5 1.25], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen});
preg.plotParamComparison(cmp, 'legend','legend', 'YLim', [-0.25 1.2], 'print','basic', 'MaxTime', 3, 'Color', {customBlue,customOrange, customGreen});

ar.triplecode.updatePlotParm({'Operand 1', 'Operand 2', 'Operand 3'}, 'A+(B+C)')

save2pdf([figuresDir 'rr_Ap_BpC_Op1Op2Op3PrevTldrfix.pdf'], gcf, 600) 












% Update the filtering (to skip the whole preprocessing, since postprocess_set() was updated after it)
% postprocess_set(calc1.raw)
% postprocess_set(calc2.raw)

%% Update y velocity and acceletation for calc1 and calc2 (Trajectory columns 16 and 17)
% [calc1] = ar.calc1.updateYVelAcc('OM_exp1', 0.02);
% [calc2] = ar.calc1.updateYVelAcc('OM_exp2', 0.02);

%% Filter trials with errors
% EndPointErrorCheck = cell2mat(arrayfun(@(t){t.Custom.EndPointPercError}, triplecode.d.all.Trials));

% Type = arrayfun(@(t){t.Custom.Type}, triplecode.d.all.Trials);

twoOpsFil = tt.util.filterDataset(twoOps, @(t) t.Custom.EndPointPercError <= .1); 

length(twoOpsFil.raw.all.Trials)/length(twoOps.raw.all.Trials)







%% Prepare manual encoding of horizontal velocity 
% ar.calc1.prepareOnsetVelocities(1, 'NoManual');








%%


















% Additions L+S
cmp = treg.compareParams(rr_c1addLS, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero');
preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions L+S')
save2pdf([figuresDir 'rr_c1addLS_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions S+L
cmp = treg.compareParams(rr_c1addSL, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero');
preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions S+L')
save2pdf([figuresDir 'rr_c1addSL_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions and Subtractions for Operational momentum with Abs operand2
cmp = treg.compareParams(rr_c1addsubm, 'ADD_SUB_Op1OperatorSignOp2Op2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2','b_Operand2','b_Operator'}, 'compare1method', 'Tzero');
preg.plotParamComparison(cmp, 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {customGrey,lightGrey,lightGrey-0.08,'black'}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2','Absolute Operand 2','Operator'}, 'Operational Momentum Effect')
save2pdf([figuresDir 'ADD_SUB_Op1OperatorSignOp2Op2prevTargetldrfix.pdf'], gcf, 600)

% % Additions and Subtractions for Operational momentum 
% preg.plotParamComparison(treg.compareParams(rr_c1addsubm, 'ADD_SUB_Op1OperatorSignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2', 'b_Operator'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {customGrey,lightGrey,'black', 'red'}, 'print','basic');
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2', 'Operator'}, 'Operational Momentum Effect')
% save2pdf([figuresDir 'rr_c1addsubm_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Control with no 0 operands

% Subtractions No 0
% preg.plotParamComparison(treg.compareParams(rr_c1subNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Subtractions with no Operands = 0')
% save2pdf([figuresDir 'rr_c1subNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions L+S No 0
% preg.plotParamComparison(treg.compareParams(rr_c1addLSNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions L+S with no Operands = 0')
% save2pdf([figuresDir 'rr_c1addLSNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions S+L No 0
% preg.plotParamComparison(treg.compareParams(rr_c1addSLNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2', 'b_ldrfix', 'b_prevtarget'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey, 'blue', 'red'}, 'print','basic');
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions S+L with no Operands = 0' )
% save2pdf([figuresDir 'rr_c1addSLNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions and Subtractions for Operational momentum No 0
% preg.plotParamComparison(treg.compareParams(rr_c1addsubmNo0, 'ADD_SUB_Op1OperatorSignOp2Op2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2', 'b_Operator' 'b_Operand2',}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {customGrey,lightGrey,'black', 'red'}, 'print','basic');
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2', 'Operator'}, 'Operational Momentum Effect')
% save2pdf([figuresDir 'rr_c1addsubmNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Control for range
% Subtractions
preg.plotParamComparison(treg.compareParams(rr_c1subRange, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% Update plot parameters and save plot
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Subtractions')
save2pdf([figuresDir 'rr_c1subRange_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions L+S
preg.plotParamComparison(treg.compareParams(rr_c1addLSRange, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions L+S')
save2pdf([figuresDir 'rr_c1addLSRange_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions S+L
preg.plotParamComparison(treg.compareParams(rr_c1addSLRange, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions S+L')
save2pdf([figuresDir 'rr_c1addSLRange_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Control for range with no O
% Subtractions
preg.plotParamComparison(treg.compareParams(rr_c1subRangeNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% Update plot parameters and save plot
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Subtractions same range no 0')
save2pdf([figuresDir 'rr_c1subRangeNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions L+S
preg.plotParamComparison(treg.compareParams(rr_c1addLSRangeNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions L+S same range no 0')
save2pdf([figuresDir 'rr_c1addLSRangeNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Additions S+L
preg.plotParamComparison(treg.compareParams(rr_c1addSLRangeNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions S+L same range no 0')
save2pdf([figuresDir 'rr_c1addSLRangeNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)


%% Comparing regression coefficients across conditions
% 
% rrAddSubCompOp1 = preg.combineRRSets({rr_c1subm, rr_c1addLS}, {'SUB', 'ADD'}, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', 'ParamName', 'Operand1');
% rrAddSubCompSignOp2 = preg.combineRRSets({rr_c1subm, rr_c1addLS}, {'SUB', 'ADD'}, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', 'ParamName', 'SignOperand2');
% rrAddLSSLCompOp1 = preg.combineRRSets({rr_c1addLS,rr_c1addSL}, {'OpMinLS', 'OpMinSL'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'OperandMin');
% rrAddLSSLCompSignOp2 = preg.combineRRSets({rr_c1addLS, rr_c1addSL}, {'OpMaxLS', 'OpMaxLS'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'OperandMax');

% Combine min and max in additions LS and SL
rr_AddLSSLCompOpMin = preg.combineRRSets({rr_c1addLSMinMax, rr_c1addSLMinMax}, {'addLSOpMin', 'addSLOpMin'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'SignOperandMin');
rr_AddLSSLCompOpMax = preg.combineRRSets({rr_c1addLSMinMax, rr_c1addSLMinMax}, {'addLSOpMax', 'addSLOpMax'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'OperandMax');

cmp=treg.compareParams(rr_AddLSSLCompOpMin, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_addLSOpMin', 'b_addSLOpMin'}, 'compare', 1:2, 'Cmp1Larger');
[cmp.times cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]
% preg.plotParamComparison(cmp, 'fill', 'box')
preg.plotParamComparison(cmp, 'YLim', [-0.25 1.25])
ar.calc1.updatePlotParm({'p <= .05', 'p <= .01'}, 'Min Operand')
text(.3,.5, 'Addition L+S', 'FontSize', 30)
text(.7,.5, 'Addition S+L', 'FontSize', 30)
save2pdf([figuresDir 'rr_c1MinComparisonAddLSAddSL.pdf'], gcf, 600)


cmp=treg.compareParams(rr_AddLSSLCompOpMax, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_addLSOpMax', 'b_addSLOpMax'}, 'compare', 1:2, 'Cmp1Larger');
[cmp.times cmp.values1 cmp.values2 cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]
% preg.plotParamComparison(cmp, 'fill', 'box')
preg.plotParamComparison(cmp, 'YLim', [-0.25 1.25])
ar.calc1.updatePlotParm({'p <= .05', 'p <= .01'}, 'Max Operand')
text(.2,.5, 'Addition L+S', 'FontSize', 30)
text(.6,.5, 'Addition S+L', 'FontSize', 30)
save2pdf([figuresDir 'rr_c1MaxComparisonAddLSAddSL.pdf'], gcf, 600)

% Combine min and max in additions and subtractions with the same range
rr_SubAddLSSLCompOpMin = preg.combineRRSets({rr_c1submMinMax,rr_c1addLSMinMax, rr_c1addSLMinMax}, {'subOpMin','addLSOpMin', 'addSLOpMin'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'SignOperandMin');
cmp=treg.compareParams(rr_SubAddLSSLCompOpMin, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMin', 'b_addLSOpMin', 'b_addSLOpMin'}, 'compare', 1:3);
[cmp.times cmp.values1 cmp.values2 cmp.values3 cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]
preg.plotParamComparison(cmp, 'YLim', [-0.25 1.25])
ar.calc1.updatePlotParm({'p <= .05', 'p <= .01', 'p <= .001'}, 'Min Operand')
hLeg = legend(gca);
set(hLeg,'visible','off')
save2pdf([figuresDir 'rr_c1MinComparisonSubAddLSAddSL.pdf'], gcf, 600)


rr_SubAddLSSLCompOpMax = preg.combineRRSets({rr_c1submMinMax,rr_c1addLSMinMax, rr_c1addSLMinMax}, {'subOpMax','addLSOpMax', 'addSLOpMax'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'OperandMax');
cmp=treg.compareParams(rr_SubAddLSSLCompOpMax, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMax', 'b_addLSOpMax', 'b_addSLOpMax'}, 'compare', 1:3);
[cmp.times cmp.values1 cmp.values2 cmp.values3 cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]
preg.plotParamComparison(cmp, 'YLim', [-0.25 1.25])
ar.calc1.updatePlotParm({'p <= .05'}, 'Max Operand')
hLeg = legend(gca);
set(hLeg,'visible','off')
save2pdf([figuresDir 'rr_c1MaxComparisonSubAddLSAddSL.pdf'], gcf, 600)

% Plot comparison min
preg.plotParamComparison(treg.compareParams(rr_SubAddLSSLCompOpMin, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMin', 'b_addLSOpMin', 'b_addSLOpMin'}, 'compare1method', 'Tzero'), 'legend','legend', 'YLim', [-0.25 1.25], 'color', {'black',customGrey, lightGrey}, 'print','basic');
preg.plotParamComparison(treg.compareParams(rr_SubAddLSSLCompOpMin, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'beta_subOpMin', 'beta_addLSOpMin', 'beta_addSLOpMin'}, 'compare1method', 'Tzero'), 'legend','legend', 'YLim', [-0.1 0.45], 'color', {'black',customGrey, lightGrey}, 'print','basic');

ar.calc1.updatePlotParm({'Subtractions', 'Additions L+S', 'Additions S+L'}, 'Min Operand')
save2pdf([figuresDir 'calc1_RegCompMin_SubAdd.pdf'], gcf, 600)

% Plot comparison max
preg.plotParamComparison(treg.compareParams(rr_SubAddLSSLCompOpMax, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMax', 'b_addLSOpMax', 'b_addSLOpMax'}, 'compare1method', 'Tzero'), 'legend','legend', 'YLim', [-0.25 1.25], 'color', {'black',customGrey, lightGrey}, 'print','basic');
preg.plotParamComparison(treg.compareParams(rr_SubAddLSSLCompOpMax, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'beta_subOpMax', 'beta_addLSOpMax', 'beta_addSLOpMax'}, 'compare1method', 'Tzero'), 'legend','legend', 'YLim', [-0.25 1.25], 'color', {'black',customGrey, lightGrey}, 'print','basic');

ar.calc1.updatePlotParm({'Subtractions', 'Additions L+S', 'Additions S+L'}, 'Max Operand')
save2pdf([figuresDir 'calc1_RegCompMax_SubAdd.pdf'], gcf, 600)

% Plot comparison min and max with stats

% Stats
% cmp=treg.compareParams(rrAddSubCompSignOp2, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_SUB', 'b_ADD'}, 'compare', 1:2, 'Cmp1Larger');
% [cmp.times cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]


%% Curvature analysis

%% Seriality analysis 2 - successive trajectories and subtraction of sucsessive trajectories
ar.calc1.display_xbyt_average('exp1', '1+')
figure(1)
save2pdf([figuresDir 'traj_9-_exp1.pdf'], gcf, 600)

ar.calc1.display_xbyt_average('exp1', '9-')
figure(1)
save2pdf([figuresDir 'traj_9-_exp1.pdf'], gcf, 600)

[timexvalue, statsReg] = ar.calc1.subtract_trajs_add_sub('exp1', '1+');
figure(1)
save2pdf([figuresDir 'subtraj_1+_exp1.pdf'], gcf, 600)
figure(2)
save2pdf([figuresDir 'subtraj_slopes_1+_exp1.pdf'], gcf, 600)
[arrayfun(@(t)t.beta(2), statsReg); arrayfun(@(t)t.p(2), statsReg)]

[timexvalue, statsReg] = ar.calc1.subtract_trajs_add_sub('exp1', '9-');
figure(1)
save2pdf([figuresDir 'subtraj_9-_exp1.pdf'], gcf, 600)
figure(2)
save2pdf([figuresDir 'subtraj_slopes_9-_exp1.pdf'], gcf, 600)
[arrayfun(@(t)t.beta(2), statsReg); arrayfun(@(t)t.p(2), statsReg)]

%% 2. Loading Exp2
calc2  = tt.loadDataset('OM_exp2/all');
% Update the filtering (to skip the whole preprocessing, since postprocess_set() was updated after it)
% postprocess_set(calc2.raw)

tt.util.filterDataset(calc2.raw, @(t)~t.Custom.isCalc, 'Save', 'Dir', 'OM_Exp2/1digit');
tt.util.filterDataset(calc2.raw, @(t)t.Custom.isCalc && t.Custom.Operator == 1, 'Save', 'Dir', 'OM_Exp2/add');
tt.util.filterDataset(calc2.raw, @(t)t.Custom.isCalc && t.Custom.Operator == -1, 'Save', 'Dir', 'OM_Exp2/sub');
tt.util.filterDataset(calc2.raw, @(t)t.Custom.isCalc,'Save', 'Dir','OM_Exp2/add_sub');

c21d = tt.loadDataset('OM_exp2/1digit');
c2add = tt.loadDataset('OM_exp2/add');
c2sub = tt.loadDataset('OM_exp2/sub');
c2addsub = tt.loadDataset('OM_exp2/add_sub');

%% Prepare manual encoding of horizontal velocity 
ar.calc1.prepareOnsetVelocities(2, 'NoManual');
% nl.vel.plotVelocityOnset(c2add.d.am); % To review the automatic encoding
% nl.vel.encodeVelocityOnsetManually({c2add.d, c2sub.d, c21d.d}, @(t)~nl.vel.trialHasOnset(t, 1, 1, 0), '/tmp/onset-add-sub_exp2.csv');

% Velocity onset
ar.calc1.prepareOnsetVelocities(2)

% Dynamic filtering
c2addNo0 = tt.util.filterDataset(c2add, @(t)t.Custom.Operand1 ~= 0 && t.Custom.Operand2 ~= 0) ; 
c2subNo0 = tt.util.filterDataset(c2sub, @(t)t.Custom.Operand1 ~= 0 && t.Custom.Operand2 ~= 0) ; 
c2addsubNo0 = tt.util.filterDataset(c2addsub, @(t)t.Custom.Operand1 ~= 0 && t.Custom.Operand2 ~= 0) ; 

%% Geneneral Performance Measures
generalPerfc2 = tt.inf.printBasicStats({c2sub,c2add});
[H,P,CI,STATS] = ttest(generalPerfc2{1}.MovementTime, generalPerfc2{2}.MovementTime);
[H,P,CI,STATS] = ttest(generalPerfc2{1}.EndpointErr, generalPerfc2{2}.EndpointErr);
[H,P,CI,STATS] = ttest(generalPerfc2{1}.EndpointBias, generalPerfc2{2}.EndpointBias);

% Averaged trajectories by result
ar.calc1.generalMeasuresPlot('OM_exp2/sub', 'OM_exp2/add',2,0);
figure(1)
save2pdf([figuresDir 'xbyt_Sub_ADD_exp2.pdf'], gcf, 600)

% Problem difficulty effect
ar.calc1.plotProblemDifficulty(2)
figure(2)
save2pdf([figuresDir 'ProblemDifficultyEffect_exp2.pdf'], gcf, 600)

% Movement time, endpoint error and endpoint bias
% ar.calc1.plotMovTimeEPerrorEPbias(2) 
% save2pdf([figuresDir 'MovTimeEPErrorEPbias_exp2.pdf'], gcf, 600)

%% Regressions

% Experiment 2
% rr_c21d = treg.runOneRegression(c21d, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);
rr_c2add = treg.runOneRegression(c2add, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);
rr_c2sub = treg.runOneRegression(c2sub, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);
rr_c2addsub = treg.runOneRegression(c2addsub, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);

rr_c2subMinMax = treg.runOneRegression(c2sub, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);
rr_c2addMinMax = treg.runOneRegression(c2add, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);

% Control for no 0 operands
rr_c2addNo0 = treg.runOneRegression(c2addNo0, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);
rr_c2subNo0 = treg.runOneRegression(c2subNo0, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);
rr_c2addsubNo0 = treg.runOneRegression(c2addsubNo0, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', regMaxTime, 'IRD', regSampRate);


%% Plotting 

% Subtractions
preg.plotParamComparison(treg.compareParams(rr_c2sub, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% Update plot parameters and save plot
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Subtractions')
save2pdf([figuresDir 'rr_c2sub_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

times = rr_c2sub.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.times;
b_prevtarget = rr_c2sub.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.b_prevtarget;
p_prevtarget = rr_c2sub.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.p_prevtarget;
b_ldrfix = rr_c2sub.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.b_ldrfix;
p_ldrfix = rr_c2sub.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.p_ldrfix;
table(times, b_prevtarget, p_prevtarget, b_ldrfix, p_ldrfix)

% Additions
preg.plotParamComparison(treg.compareParams(rr_c2add, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% Update plot parameters and save plot
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions L+S')
save2pdf([figuresDir 'rr_c2addLS_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

times = rr_c2add.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.times;
b_prevtarget = rr_c2add.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.b_prevtarget;
p_prevtarget = rr_c2add.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.p_prevtarget;
b_ldrfix = rr_c2add.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.b_ldrfix;
p_ldrfix = rr_c2add.avg.ADD_OR_SUB_Op1SignOp2prevTargetldrfix.p_ldrfix;
table(times, b_prevtarget, p_prevtarget, b_ldrfix, p_ldrfix)

% Additions and Subtractions for Operational momentum 
preg.plotParamComparison(treg.compareParams(rr_c2addsub, 'ADD_SUB_Op1OperatorSignOp2Op2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2','b_Operand2','b_Operator'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {customGrey,lightGrey,lightGrey-0.08,'black'}, 'print','basic');
% Update plot parameters and save plot
ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2', 'Absolute Operand 2','Operator'}, 'Operational Momentum Effect')
save2pdf([figuresDir 'rr_c2addsub_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)

% Control for no 0 operands
% Subtractions No 0
% preg.plotParamComparison(treg.compareParams(rr_c2subNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% % Update plot parameters and save plot
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Subtractions with no Operands = 0')
% save2pdf([figuresDir 'rr_c2subNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)
% 
% % Additions no 0
% preg.plotParamComparison(treg.compareParams(rr_c2addNo0, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
% % Update plot parameters and save plot
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2'}, 'Additions L+S with no Operands = 0')
% save2pdf([figuresDir 'rr_c2addNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)
% 
% % Additions and Subtractions for Operational momentum No 0
% preg.plotParamComparison(treg.compareParams(rr_c2addsubNo0, 'ADD_SUB_Op1OperatorSignOp2Op2prevTargetldrfix', {'b_Operand1', 'b_SignOperand2', 'b_Operator' 'b_Operand2',}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {customGrey,lightGrey,'black', 'red'}, 'print','basic');
% % Update plot parameters and save plot
% ar.calc1.updatePlotParm({'Operand 1', 'Signed Operand 2', 'Operator'}, 'Operational Momentum Effect with no Operands = 0')
% save2pdf([figuresDir 'rr_c2addsubNo0_Op1SignOp2prevTargetldrfix.pdf'], gcf, 600)
% 

%% Comparing regression coefficients across conditions
%
% Combine min and max in additions and subtractions
rrcalc2SubAddCompOpMin = preg.combineRRSets({rr_c2subMinMax,rr_c2addMinMax}, {'subOpMin','addOpMin'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'SignOperandMin');
rrcalc2SubAddCompOpMax = preg.combineRRSets({rr_c2subMinMax,rr_c2addMinMax}, {'subOpMax','addOpMax'}, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', 'ParamName', 'OperandMax');

% Min
cmp=treg.compareParams(rrcalc2SubAddCompOpMin, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMin', 'b_addOpMin'}, 'compare', 1:2, 'Cmp1Larger');
[cmp.times cmp.values1 cmp.values2 cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]
preg.plotParamComparison(cmp, 'YLim', [-0.25 1.25])
ar.calc1.updatePlotParm({'p <= .05', 'p <= .01', 'p <= .001'}, 'Min Operand')
hLeg = legend(gca);
set(hLeg,'visible','off')
% text(.2,.5, 'Addition L+S', 'FontSize', 30)
% text(.6,.5, 'Addition S+L', 'FontSize', 30)
save2pdf([figuresDir 'rr_c2MinComparisonAdd.pdf'], gcf, 600)

% Max
cmp=treg.compareParams(rrcalc2SubAddCompOpMax, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMax', 'b_addOpMax'}, 'compare', 1:2, 'Cmp1Larger');
[cmp.times cmp.values1 cmp.values2 cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]
preg.plotParamComparison(cmp, 'YLim', [-0.25 1.25])
ar.calc1.updatePlotParm({'p <= .05'}, 'Max Operand')
hLeg = legend(gca);
set(hLeg,'visible','off')
% text(.2,.5, 'Addition L+S', 'FontSize', 30)
% text(.6,.5, 'Addition S+L', 'FontSize', 30)
save2pdf([figuresDir 'rr_c2MaxComparisonAdd.pdf'], gcf, 600)


% Plot comparison min
preg.plotParamComparison(treg.compareParams(rrcalc2SubAddCompOpMin, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMin', 'b_addOpMin'}, 'compare1method', 'Tzero'), 'legend','legend', 'YLim', [-0.25 1.25], 'color', {'black',customGrey, lightGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Subtractions', 'Additions'}, 'Min Operand')
save2pdf([figuresDir 'rr_c1MinComparisonSubAdd.pdf'], gcf, 600)

% Plot comparison max
preg.plotParamComparison(treg.compareParams(rrcalc2SubAddCompOpMax, 'ADD_OR_SUB_OpMinOpMaxTargetldrfix', {'b_subOpMax', 'b_addOpMax'}, 'compare1method', 'Tzero'), 'legend','legend', 'YLim', [-0.25 1.25], 'color', {'black',customGrey, lightGrey}, 'print','basic');
ar.calc1.updatePlotParm({'Subtractions', 'Additions'}, 'Max Operand')
save2pdf([figuresDir 'rr_c1MaxComparisonSubAdd.pdf'], gcf, 600)


%% END



























%%



% Stats
% cmp=treg.compareParams(rrAddSubCompSignOp2, 'ADD_OR_SUB_Op1SignOp2prevTargetldrfix', {'b_SUB', 'b_ADD'}, 'compare', 1:2, 'Cmp1Larger');
% [cmp.times cmp.comparePair.pParam cmp.comparePair.fParam cmp.comparePair.dfParam]

%%

%% Single digit preceeded by addition or subtractions
preg.plotParamComparison(treg.compareParams(rr_c1sub, 'ADD_OR_SUB_Op1SignOp2', 'allb', 'compare1method', 'Tzero'),'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic')
preg.plotParamComparison(treg.compareParams(rr_c1subm, 'ADD_OR_SUB_Op1SignOp2', 'allb', 'compare1method', 'Tzero'),'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic')
rrPlugOp1 = preg.combineRRSets({rr_c1sub, rr_c1subm}, {'c1sub', 'c1subm'}, 'ADD_OR_SUB_Op1SignOp2', 'ParamName', 'Operand1');
rrPlugSignOp2 = preg.combineRRSets({rr_c1sub, rr_c1subm}, {'c1sub', 'c1subm'}, 'ADD_OR_SUB_Op1SignOp2', 'ParamName','SignOperand2');
preg.plotParamComparison(treg.compareParams(rrPlugOp1, 'ADD_OR_SUB_Op1SignOp2', {'b_c1sub', 'b_c1subm'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic');
preg.plotParamComparison(treg.compareParams(rrPlugSignOp2, 'ADD_OR_SUB_Op1SignOp2', {'b_c1sub', 'b_c1subm'}, 'compare1method', 'Tzero'), 'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {customGrey, 'grey', 'black'}, 'print','basic');
ar.calc1.updatePlotParm({'c1sub', 'c1sub'})




save2pdf([figuresDir 'test_0+D_D_effect_of_target.pdf'], gcf, 600)

%% Single digit preceeded by subtraction




%%

preg.plotParamComparison(treg.compareParams(rr_c1add0D, 'SD_targetRefPoint3', 'allb', 'compare1method', 'Tzero'),...
    'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic')

preg.plotParamComparison(treg.compareParams(rr_c11d, 'SD_targetRefPoint3', 'allb', 'compare1method', 'Tzero'),...
    'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey}, 'print','basic')


rrpnDec = preg.combineRRSets({rr_c1add0D, rr_c1addD0}, {'0D', 'D0'}, 'SD_targetRefPoint3', 'ParamName', 'Target');


preg.plotParamComparison(treg.compareParams(rrpnDec, 'SD_targetRefPoint3', {'b_0D', 'b_D0'}), 'LineStyle', {'-', '-'}, 'Color', {'blue', 'red'}, 'Marker', {'o', 'd'}, 'MarkerSize', 20, 'YLim', [-.2 1.45], 'Legend', 'none','MaxTime', 1);
save2pdf([figuresDir 'test_0+D_D_effect_of_target.pdf'], gcf, 600)



preg.plotParamComparison(treg.compareParams(rr_c1addLS, 'ADD_OR_SUB_Op1SignOp2prevTarget', {'b_SignOperand2'}), 'LineStyle', {'-', '-'}, 'Color', {'blue', 'red'}, 'Marker', {'o', 'd'}, 'MarkerSize', 20, 'YLim', [-.2 1.45], 'Legend', 'None', 'MaxTime', 1);


%%




%% Prepare manual encoding of horizontal velocity 
ar.calc1.prepareOnsetVelocities(1, 'NoManual');
nl.vel.plotVelocityOnset(c1add.d.as); % To review the automatic encoding
nl.vel.encodeVelocityOnsetManually({c1add.d, c1sub.d, c11d.d}, @(t)~nl.vel.trialHasOnset(t, 1, 1, 0), '/tmp/onset-add-sub.csv');

% Velocity onset
ar.calc1.prepareOnsetVelocities(1)


% tt.util.filterDataset(Data_OM_exp1.raw, @(t)t.Custom.isCalc && t.Custom.Operator == 1 && t.Custom.Operand1 > t.Custom.Operand2, 'OM_Exp1/addLS_matched_op');


% Curvature analysis









rr_c1sub = treg.runOneRegression(c1sub, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', 1.31, 'IRD', regSampRate);

preg.plotParamComparison(treg.compareParams(rr_c1sub, 'ADD_OR_SUB_Op1SignOp2prevTarget', 'allb', 'compare1method', 'Tzero'),...
    'legend','legend', 'ci', 0.95, 'YLim', [-0.25 1.25], 'color', {'black',customGrey, 'blue', 'blue'}, 'print','basic')


preg.plotParamComparison(treg.compareParams(rr_c1sub, 'ADD_OR_SUB_Op1SignOp2prevTarget', {'Operand1','SignOperand2','prevtarget'}),...
    'LineStyle', {'-', '--', '-'}, 'Color', {'blue', 'blue', 'red'}, 'Marker', {'o', 'd', 'o'},...
    'MarkerSize', 20, 'YLim', [-.2 1.45], 'Legend', 'None', 'MaxTime', 1);



% Correct legend names
hLeg = legend(gca);
hLeg.String = {'Operand 1', 'Signed operand 2'};
set(gca, 'XTick', 0:0.25:1.5);
set(gca, 'XTickLabel',  0:0.25:1.5);
set(gcf, 'Position', [0 0 900*1.3 600*1.3])
textfig = findobj(gcf);
textfig = findall(textfig,'Type','text');
delete(textfig);
ylabel('b','rot',360)
xlabel('Time (sec.)')
ylabh = get(gca,'YLabel');
xlabh = get(gca,'XLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.01 0 0])
% title('Exp1 ADD L+S Op1SignOp2','FontSize', FontSizeTitle)




% rr_AddSLsmall_Exp1 = treg.runOneRegression(out_addSLsmall, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', 1.31, 'IRD', regSampRate);
% rr_AddSLlarge_Exp1 = treg.runOneRegression(out_addSLlarge, 'ep', 'RemoveOutliers', 'ep', 'AllTrials', 'MaxTime', 1.31, 'IRD', regSampRate);

% 





common.alignreg.tryAlignOneRR({rrpnf00.avg.decadeUnitLDRPtrg, rrpnf10.avg.decadeUnitLDRPtrg}, 'b_decades', [.35 .65], -.2, .2, 'CondNames', {'0,0', '100,0'}, 'print');
