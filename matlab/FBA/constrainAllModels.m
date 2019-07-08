function [unconstrainedF, constrainedF] = constrainAllModels(modelsDir)

    numCellLines = 96;
    unconstrainedF = zeros(96,1);
    constrainedF = zeros(96,1);
    
    for i = 1:numCellLines
        modelFile = sprintf(modelsDir + "model_P6975_%d.mat", (i + 100));
        
        choModel = load(modelFile);
        choModel = choModel.results.tissueModel;
        constrainedModel = setKallehaugeBounds(choModel);
        
        unconstrainedFBA = optimizeCbModel(choModel);
        constrainedFBA = optimizeCbModel(constrainedModel);
    
        unconstrainedF(i) = unconstrainedFBA.f;
        constrainedF(i) = constrainedFBA.f;
    end
end
