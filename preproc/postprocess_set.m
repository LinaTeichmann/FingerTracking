%
% This script is invoked by the preprocessSet() function
%
function raw = postprocess_set(raw)
    % Mark decade and unit delays
    subjIDs = tt.inf.listInitials(raw);
    for i = 1:length(subjIDs)
        expData = raw.(subjIDs{i});
        prevOperand1 = [];
        prevOperand2 = [];
        prevOperator = [];
        prevResTarget = [];        
        for trial = expData.Trials
            trial.Custom.prevOperand1 = prevOperand1;
            trial.Custom.prevOperand2 = prevOperand2;
            trial.Custom.prevOperator = prevOperator;
            trial.Custom.prevResTarget = prevResTarget;
            trial.Custom.ResTarget = str2num(trial.Custom.PresentedTarget);
            Operands = strsplit(trial.Custom.PresentedTarget, {'-', '+'});            
            trial.Custom.Operand1 = str2num(Operands{1});
            trial.Custom.Operand2 = str2num(Operands{2});
            trial.Custom.OperandMin = min([trial.Custom.Operand1, trial.Custom.Operand2]);
            trial.Custom.OperandMax = max([trial.Custom.Operand1, trial.Custom.Operand2]);
            trial.Custom.DistOps = abs(trial.Custom.Operand1 - trial.Custom.Operand2); 
            trial.Custom.RatioOps = trial.Custom.OperandMin/trial.Custom.OperandMax;
            trial.Custom.EndPointPercError = 1-(min(trial.EndPoint, trial.Target)/max(trial.EndPoint, trial.Target));       
            trial.Custom.Product = trial.Custom.Operand1*trial.Custom.Operand2; 

            OperatorIdx = strfind(trial.Custom.PresentedTarget, '-');
                if isempty(OperatorIdx) == 1
                    trial.Custom.Operator = 1;
                else
                    trial.Custom.Operator = -1; % subtraction
                end
                
                if trial.Custom.Operator == 1 && trial.Custom.Operand1 > trial.Custom.Operand2
                   trial.Custom.Type = 'L+S';
                elseif trial.Custom.Operator == 1 && trial.Custom.Operand1 < trial.Custom.Operand2
                   trial.Custom.Type = 'S+L';
                else
                   trial.Custom.Type = 'L-S';
                end
                
            prevOperand1 = trial.Custom.Operand1;
            prevOperand2 = trial.Custom.Operand2;
            prevOperator = trial.Custom.Operator;
            prevResTarget = trial.Custom.ResTarget;
        end
    end
end