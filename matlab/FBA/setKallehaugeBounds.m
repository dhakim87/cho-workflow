function [alteredModel] = setKallehaugeBounds(choModel)
    constraints = load('experimental_constraints.mat');
    alteredModel = choModel;
    
    bounds = constraints.kallehauge_bounds;
    reactions = constraints.kallehauge_rxns;
    manuallyUpdatedReactions = [...
        cellstr('EX_ala__L_e'), cellstr('EX_arg__L_e'), cellstr('EX_asn__L_e'), ...
        cellstr('EX_asp__L_e'), cellstr('EX_cys__L_e'), cellstr('EX_gln__L_e'), ...
        cellstr('EX_glu__L_e'), cellstr('EX_gly_e'),    cellstr('EX_his__L_e'), ...
        cellstr('EX_ile__L_e'), cellstr('EX_leu__L_e'), cellstr('EX_lys__L_e'),...
        cellstr('EX_met__L_e'), cellstr('EX_phe__L_e'), cellstr('EX_pro__L_e'),...
        cellstr('EX_ser__L_e'), cellstr('EX_thr__L_e'), cellstr('EX_trp__L_e'),...
        cellstr('EX_tyr__L_e'), cellstr('EX_val__L_e'), cellstr('EX_glc__D_e'), ...
        cellstr('EX_lac__L_e')];
    
    numReactionsFound = 0;
    numReactionsNotFound = 0;
    for kallehaugeReactionIndex = 1:numel(reactions)
%         kalleReaction = reactions{kallehaugeReactionIndex};
        updatedString = manuallyUpdatedReactions{kallehaugeReactionIndex};
        foundIndex = 0;
        numFound = 0;
%         fprintf("Reaction:         %s\n", kalleReaction);
%         fprintf("Manually Updated: %s\n", updatedString);
        for choReactionIndex = 1:numel(choModel.rxns)
            choReaction = choModel.rxns{choReactionIndex};
            if (strcmp(updatedString, choReaction))
%                 fprintf('\tCHO Reaction: %s\n', choReaction);
                foundIndex = choReactionIndex;
                numFound = numFound + 1;
            end
        end
%         if (~foundIndex)
%             fprintf('Could not find: %s\n', updatedString);
%         end
        if (numFound > 1)
            fprintf("WARNING:Found more than one answer!\n");
        end
        if (numFound == 1)
            numReactionsFound = numReactionsFound + 1;
        else
            numReactionsNotFound = numReactionsNotFound + 1;
        end
%         fprintf("Found Reaction Index: %d\n", foundIndex);
    end
    
    fprintf("Num Reactions Found    : %d\n", numReactionsFound);
    fprintf("Num Reactions Not Found: %d\n", numReactionsNotFound); 
    
    for kallehaugeReactionIndex = 1:numel(reactions)
        choReactionName = manuallyUpdatedReactions{kallehaugeReactionIndex};
        lb = bounds(kallehaugeReactionIndex, 1);
        ub = bounds(kallehaugeReactionIndex, 2);
        
        if (lb > ub)
           tmp = lb;
           lb = ub;
           ub = tmp;
        end
        
        alteredModel = changeRxnBounds(alteredModel, choReactionName, lb, 'l');
        alteredModel = changeRxnBounds(alteredModel, choReactionName, ub, 'u');
    end
    
%     if isfield(alteredModel, 'csense')
%         [~, cCols] = size(alteredModel.csense);
%         if (cCols == 1)
%             alteredModel.csense = transpose(alteredModel.csense);
%         end
%     end
end