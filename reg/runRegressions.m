function rr = runRegressions(allExpData, varargin)
%rr = runRegressions(allExpData) - run regressions for all subjects.
% 
% Copy this function to your own folder and customize it.
% Note all places with a "CUSTOM" comment

    [outFN, regArgs] = parseArgs(varargin);

    if isfield(allExpData, 'raw') 
        allED = allExpData.d; 
    else
        allED = allExpData;
    end
    edArray = tt.util.structToArray(allED);
    
    rr = tt.reg.createEmptyRR('NL', allED.general.setName, edArray(1).MaxTarget);
    
    for expData = edArray
        rr.(expData.SubjectInitials) = runForOneSubj(expData, regArgs);
    end
    
    rr.avg = tt.reg.averageRegressionResults(rr);
    
    if ~isempty(outFN)
        outDirName = [TrajTrackerDataPath '/arit_express/results/regressions/', outFN];
        save(outDirName, 'rr');
    end
    
    %------------------------------------------------------------------------
    function result = runForOneSubj(expData, regArgs)
        
        result = tt.reg.createEmptyRROneSubj(expData.SubjectInitials, expData.SubjectName, expData.SamplingRate);
        result.MaxMovementTime = max(arrayfun(@(t)t.MovementTime, expData.Trials));
        
        %% CUSTOM: change all lines below according to your regressions
%            result.Op123 = tt.reg.regress(expData, 'reg', 'ep', {'Operand1', 'Operand2', 'Operand3','prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
%            result.OpsProdPar = tt.reg.regress(expData, 'reg', 'ep', {'Operand1', 'Operand2', 'Operand3', 'ProdParent', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
%            result.OutSumProdPar = tt.reg.regress(expData, 'reg', 'ep', {'OpOutParent', 'SumParent', 'ProdParent', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
%            result.OutMaxMinPar = tt.reg.regress(expData, 'reg', 'ep', {'OpOutParent', 'OpMaxParent', 'OpMinParent', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
%            result.OutPar = tt.reg.regress(expData, 'reg', 'ep', {'OpOutParent', 'ResParent', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
%            result.OutTarget = tt.reg.regress(expData, 'reg', 'ep', {'OpOutParent', 'ResTarget', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);

                    
           switch(allED.general.CondName)
               case {'(AxB)+C', 'A+(BxC)', '(AxB)-C'}
                   %                  result.OutSumProdParMinMaxPar = tt.reg.regress(expData, 'reg', 'ep', {'OpOutParent', 'SumParent', 'ProdParent', 'OpMaxParent', 'OpMinParent', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
               case {'(A+B)+C', 'A+(B+C)','(A+B)-C'}
                   %                  result.OutProdParMinMaxPar = tt.reg.regress(expData, 'reg', 'ep', {'OpOutParent', 'ProdParent', 'OpMaxParent', 'OpMinParent', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
               case {'add_0_20', 'mult_0_20'}
                   result.Op12 = tt.reg.regress(expData, 'reg', 'ep', {'Operand1', 'Operand2','prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
                   result.MaxMin = tt.reg.regress(expData, 'reg', 'ep', {'OperandMin', 'OperandMax', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
                   result.Target = tt.reg.regress(expData, 'reg', 'ep', {'ResTarget', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
               case {'(A+B)+/-C'}
%                    result.Op12sign3abs3 = tt.reg.regress(expData, 'reg', 'ep', {'Operand1', 'Operand2', 'Operand3', 'AbsOperand3', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
                   result.Op12sign3abs3operator2 = tt.reg.regress(expData, 'reg', 'ep', {'Operand1', 'Operand2', 'Operand3', 'AbsOperand3', 'Operator2', 'prevtarget', 'ldrfix'}, 'TPdep','FMeasureFunc', @getMyTrialMeasures,regArgs);
               otherwise
                   error('No filtered data with this name')
           end
        
        
    end

    %-------------------------------------------
    function [outFN, regArgs] = parseArgs(args)

        outFN = 'default_filename.mat'; %CUSTOM
        regArgs = {};
        
        args = stripArgs(args);
        while ~isempty(args)
            switch(lower(args{1}))
                case 'nosave'
                    outFN = '';

                case 'saveas'
                    outFN = args{2};
                    args = args(2:end);

                case 'regargs'
                    regArgs = args{2};
                    args = args(2:end);
                    
                otherwise
                    error('Unsupported argument "%s"!', args{1});
            end
            args = stripArgs(args(2:end));
        end

    end

end

