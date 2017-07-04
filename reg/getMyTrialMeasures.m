function [measures, outMeasureNames, measureDescs] = getMyTrialMeasures(expData, trials, measureNames)
% [measures, outMeasureNames, measureDescs] = getMyTrialMeasures(expData, trials, measureNames) -
% Get values of regression measures that have one value per trial.
%
% trials - a column vector of trials (from the given expData)
% measureNames - cell array of measure names (see function for details)
%
% Return values -
% measures: matrix with cols = measures, rows = trials
% varNames: symbol name of each measure (cell array)
% measureDescs: text description of each measure (cell array)

    if (size(trials, 2) > 1)
        error('"trials" should be a column vector, but it has %d columns!', size(trials, 2));
    end
    
    measures = [];
    outMeasureNames = measureNames;
    measureDescs = cell(1, length(measureNames));

    for iMeasure = 1:length(measureNames)
        
        currMeasure = []; %#ok<NASGU>
        
        [measureName, measureArgs] = tt.reg.internal.parseMeasureName(measureNames{iMeasure}); %#ok<ASGLU>
        currMeasureDesc = ''; %#ok<NASGU>
        
        switch(lower(measureName))
            
            %-------- Stimulus ------------
            
            case 'operand1'
                currMeasure = arrayfun(@(t)t.Custom.Operand1, trials); %CUSTOM: replace this with a custom code that sets currMeasure to be a column vector
                currMeasureDesc = 'A'; %CUSTOM: replace this with a custom code that sets measure short description (for figures)
              
            case 'operand2'
                currMeasure = arrayfun(@(t)t.Custom.Operand2, trials); 
                currMeasureDesc = 'B'; 
                
            case 'operand3'
                currMeasure = arrayfun(@(t)t.Custom.SignOperand3, trials); 
                currMeasureDesc = 'C';  
                
            case 'absoperand3'
                currMeasure = arrayfun(@(t)t.Custom.Operand3, trials); 
                currMeasureDesc = 'C'; 
                
            case 'product12'
                currMeasure = arrayfun(@(t)t.Custom.Product12, trials); %CUSTOM: replace this with a custom code that sets currMeasure to be a column vector
                currMeasureDesc = 'AxB'; %CUSTOM: replace this with a custom code that sets measure short description (for figures)
              
            case 'product13'
                currMeasure = arrayfun(@(t)t.Custom.Product13, trials); 
                currMeasureDesc = 'AxC'; 
                
            case 'product23'
                currMeasure = arrayfun(@(t)t.Custom.Product23, trials); 
                currMeasureDesc = 'BxC'; 
                
            case 'operator1'
                currMeasure = arrayfun(@(t)t.Custom.Operator1, trials); 
                currMeasureDesc = 'Operator 1'; 
                
            case 'operator2'
                currMeasure = arrayfun(@(t)t.Custom.Operator2, trials); 
                currMeasureDesc = 'Operator 2';        
                
            case 'resparent'
                currMeasure = arrayfun(@(t)t.Custom.ResParent, trials); 
                currMeasureDesc = 'result()';
                
            case 'sumparent'
                currMeasure = arrayfun(@(t)t.Custom.OpMinParent + t.Custom.OpMaxParent, trials); 
                currMeasureDesc = 'sum()';
                
            case 'prodparent'
                currMeasure = arrayfun(@(t)t.Custom.OpMinParent * t.Custom.OpMaxParent, trials); 
                currMeasureDesc = 'prod()';
                
            case 'opoutparent'
                currMeasure = arrayfun(@(t)t.Custom.SignOpOutParent, trials); 
                currMeasureDesc = 'out()'; 
                
            case 'opminparent'
                currMeasure = arrayfun(@(t)t.Custom.OpMinParent, trials); 
                currMeasureDesc = 'min()';
                
            case 'opmaxparent'
                currMeasure = arrayfun(@(t)t.Custom.OpMaxParent, trials); 
                currMeasureDesc = 'max()';
                
            case 'operandmin'
                currMeasure = arrayfun(@(t)t.Custom.OperandMin, trials); 
                currMeasureDesc = 'Min Op';
                
            case 'operandmax'
                currMeasure = arrayfun(@(t)t.Custom.OperandMax, trials); 
                currMeasureDesc = 'Max Op';
                
            case 'restarget'
                currMeasure = arrayfun(@(t)t.Target, trials); 
%                 currMeasure = arrayfun(@(t)t.Custom.ResTarget, trials); 
                currMeasureDesc = 'Result';
                
                
                
            otherwise
                %-- Call default function
                [currMeasure, vn, currMeasureDesc] = tt.reg.getTrialMeasures(expData, trials, measureNames(iMeasure));
                outMeasureNames{iMeasure} = vn{1};
        end
        
        if isempty(currMeasure)
            error('Huh? the measure was not calcualted!');
        end
        
        
        measures = [measures currMeasure]; %#ok<*AGROW>
        measureDescs{iMeasure} = iif(isempty(currMeasureDesc), measureName, currMeasureDesc);
        
    end
    
end
