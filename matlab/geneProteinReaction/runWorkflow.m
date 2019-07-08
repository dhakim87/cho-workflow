function [expressionRxns, parsedGPR, gene_used] = runWorkflow()
    %Step 0, Load Data
    %TODO FIXME HACK:  Pass in which cell line to get expression data for,
    %or allow us to loop over them in further steps.
    [expressionData, choModel] = loadData(); 

    %Step 1, GPR Mapping
    minSum = false; %Decide whether we like min or sum better for our models
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(choModel, expressionData, minSum); 

    %Step 2, Thresholding
    
    %Step 3, Run mCADRE

    %Step 4, Validation

    %Step 5, FBA
end