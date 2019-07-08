function [expressionData, choModel] = loadData()

    tpmFile = "sampleTPM.csv";
    choModelFile = char("iCHOv1.mat");
    
    opts = detectImportOptions(tpmFile);
    varNames = opts.VariableNames ; % variable names
    varTypes = opts.VariableTypes ;
    varTypes(1) = {'char'};
    opts = setvartype(opts,varNames,varTypes);
    tpmData = readtable(tpmFile, opts);    
    
    choModel = readCbModel(choModelFile);
        
    expressionData.gene = tpmData{:,1};
    expressionData.value = tpmData{:,3};
end