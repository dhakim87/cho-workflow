function [fluxConsistentModel, fluxConsistentUbiquityScores] = createFluxConsistentModel(choModel, ubiquityScores)

%Why does verifyModel complain about this model?  Was it constructed
%incorrectly?  Is this a tool version mismatch?
choModel2 = choModel;

[numRows, ~] = size(choModel.grRules);
choModel2.rules = cell(numRows, 1);

%What dumbness.  You must replace grRules with rules.  This is done by
%creating a map from gene name to index, then parsing the grRules strings
%and rewriting them.  

% genes: 101, 102, 103, 104, 105
% grRule: 101 and 103 or 104  ->  rule: x(1) & x(3) | x(4)
[numGenes, ~] = size(choModel.genes);
geneMap = containers.Map(choModel.genes, 1:numGenes);

for i = 1:numRows
    grRule = choModel.grRules{i};
    if (numel(grRule) == 0)
        choModel2.rules{i} = '';
        continue;
    end
    
    tokens = strsplit(grRule);
    finalTokens = cell(numel(tokens), 1);
    
%     disp("Evaluating Rule!");
%     disp(grRule);
%     disp(numel(tokens));
    for tokenIndex = 1:numel(tokens)
        token = tokens{tokenIndex};
        if (numel(token) == 0)
            continue;
        end
        firstNonParen = numel(token);
        lastNonParen = 0;
%         disp("Evaluating Token!");
%         disp(token);
        for j = 1:numel(token)
            c = token(j);
            if (~(c == '(' || c == ')'))
                firstNonParen = min(firstNonParen, j);
                lastNonParen = max(lastNonParen, j);
            end
        end
        parenPrefix = token(1:firstNonParen-1);
        noParensToken = token(firstNonParen:lastNonParen);
        parenSuffix = token(lastNonParen+1:numel(token));

        replacement = '';
        if (strcmp(noParensToken, 'or'))
            replacement = '|';
        elseif (strcmp(noParensToken, 'and'))
            replacement = '&';
        elseif (isKey(geneMap, noParensToken))
            replacement = sprintf('x(%d)', geneMap(noParensToken));
        else
            error("Gene not found: " + noParensToken);
        end
        
        finalString = strcat(parenPrefix, replacement, parenSuffix);
        finalTokens{tokenIndex} = finalString;
    end
    
    replacementString = strjoin(finalTokens);
    choModel2.rules{i} = replacementString;
end


choModel2.csense = transpose(choModel.csense);

choModel3 = choModel2;

[fluxConsistentMetBool, fluxConsistentRxnBool, ~, ~, choModel4] = findFluxConsistentSubset(choModel3);

%The definition of flux consistent metabolites is unclear.  It appears that
%removing rows of S corresponding to flux inconsistent metabolites corrupts
%the matrix by allowing reactions to use unlimited amounts of these
%inconsistent metabolites.  

% fluxConsistentModel.mets = choModel2.mets(fluxConsistentMetBool);
fluxConsistentModel.mets = choModel2.mets(:);
% fluxConsistentModel.metNames = choModel2.metNames(fluxConsistentMetBool);
fluxConsistentModel.metNames = choModel2.metNames(:);
% fluxConsistentModel.metFormulas = choModel2.metFormulas(fluxConsistentMetBool);
fluxConsistentModel.metFormulas = choModel2.metFormulas(:);
% fluxConsistentModel.metCharge = choModel2.metCharge(fluxConsistentMetBool);
fluxConsistentModel.metCharge = choModel2.metCharge(:);

fluxConsistentModel.genes = choModel2.genes;
fluxConsistentModel.rxnGeneMat = choModel2.rxnGeneMat(fluxConsistentRxnBool, :);
fluxConsistentModel.grRules = choModel2.grRules(fluxConsistentRxnBool);
fluxConsistentModel.rxns = choModel2.rxns(fluxConsistentRxnBool);
fluxConsistentModel.rxnNames = choModel2.rxnNames(fluxConsistentRxnBool);
fluxConsistentModel.subSystems = choModel2.subSystems(fluxConsistentRxnBool);
% fluxConsistentModel.S = choModel2.S(fluxConsistentMetBool, fluxConsistentRxnBool);
fluxConsistentModel.S = choModel2.S(:, fluxConsistentRxnBool);

fluxConsistentModel.lb = choModel2.lb(fluxConsistentRxnBool);
fluxConsistentModel.ub = choModel2.ub(fluxConsistentRxnBool);
% fluxConsistentModel.b = choModel2.b(fluxConsistentMetBool);
fluxConsistentModel.b = choModel2.b(:);

fluxConsistentModel.c = choModel2.c(fluxConsistentRxnBool);
fluxConsistentModel.rev = choModel2.rev(fluxConsistentRxnBool);
fluxConsistentModel.rules = choModel2.rules(fluxConsistentRxnBool);

fluxConsistentUbiquityScores = ubiquityScores(fluxConsistentRxnBool, :);

end