function [expData, choModel] = loadForStanDep()
    
    tpmFile = "sampleTPM.csv";
    choModelFile = char("iCHOv1.mat");

    opts = detectImportOptions(tpmFile);
    varNames = opts.VariableNames ; % variable names
    varTypes = opts.VariableTypes ;
    varTypes(1) = {'char'};
    opts = setvartype(opts,varNames,varTypes);
    tpmData = readtable(tpmFile, opts);    
    
    choModel = readCbModel(choModelFile);
        
    %Must create inputs for the StanDep getModel function.
    % % expressionData: a structure containing the entire expression data
        % % gene: list of gene names
        % % value: a matrix whose each column
            % % contains the expression of genes across a tissue/context and 
            % % each row contains expression of a gene across
            % % all tissues/contexts.
        % % genesymbols (optional): any alternative gene names that the
            % % user may want to keep track of.
        % % Tissue: names of conditions
    expData.gene = tpmData{1:end, 1};
    expData.Tissue = tpmData.Properties.VariableNames(3:end);
    expData.valuebyTissue = tpmData{1:end, 3:end};
end

