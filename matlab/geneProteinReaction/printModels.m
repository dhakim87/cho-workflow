function [X] = printModels(baseDir)

if (nargin == 0)
    baseDir = "modelResults";
    disp("No base dir passed, using: " + baseDir);
end 
%Print the models in the modelResults directory as csv lists of kept
%reactions.  

proteinYields = readtable("./proteinYields.csv");
for i = 1:95
    filename = sprintf(baseDir + "/model_P6975_%d.mat", i+100);
    disp(filename);
    stuff = load(filename);
    core = stuff.results.coreRxn;
    kept = stuff.results.tissueModel.rxns;
    [numCore, ~] = size(core);
    
    X(i,:) = [sprintf("P6975_%d",i+100), numCore, proteinYields{i,2}, proteinYields{i,3}];
    
    fname = sprintf(baseDir + "/P6975_%d_core.txt",i+100);
    outFile = fopen(fname, 'w');
    fprintf(outFile, "%f, %f\n", proteinYields{i,2}, proteinYields{i,3});
    for j = 1:length(core)
       fprintf(outFile, "%s\n", core{j});
    end
    fclose(outFile);
    
    fname = sprintf(baseDir + "/P6975_%d_kept.txt",i+100);
    outFile = fopen(fname, 'w');
    fprintf(outFile, "%f, %f\n", proteinYields{i,2}, proteinYields{i,3});
    for j = 1:length(kept)
       fprintf(outFile, "%s\n", kept{j});
    end
    fclose(outFile);
    
end

writematrix(X);

