function [biomass, reactionFlux, secretionYield] = evalModelFlux(baseDir, reaction)

if (nargin == 0)
    baseDir = "modelResults";
end 

disp("Searching For Models In: " + baseDir);

%Print the models in the modelResults directory as csv lists of kept
%reactions.  

reactionFlux = zeros(95,1);
biomass = zeros(95,1);
secretionYield = readtable("./proteinYields.csv");
numDeadLines = 0;
for i = 1:95
    filename = sprintf(baseDir + "/model_P6975_%d.mat", i+100);
    disp(filename);
    stuff = load(filename);
    model = stuff.results.tissueModel;
    
    reactionIndices = find(contains(model.rxns, reaction));
    fba = optimizeCbModel(model);
    if (fba.f == 0)
        fprintf("%s is dead :(\n", filename);
        numDeadLines = numDeadLines + 1;
    end
    
    if numel(reactionIndices) == 1
        reactionFlux(i) = fba.x(reactionIndices);
    else
        reactionFlux(i) = -9999;
    end
    biomass(i) = fba.f;
end

fprintf("Num Dead Lines %d\n", numDeadLines);
