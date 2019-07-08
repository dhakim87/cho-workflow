function generateCellLineModels(cellLineStartIndex, cellLineEndIndex)
    %Wrapper script to run everything on the server

    %Oh hell, the matlab path on this server is all kinds of screwed up.
    %So I have to go in and remove everybody else's paths before I can do anything.
    %First we have to fix the path, then we can actually work on doing things.
    %Need to ask Nate/others about fixing this more permanently.

    lines = strsplit(path,':');
    omarLines = startsWith(lines, "/home/omar");
    omarLines = lines(omarLines);
    [~, numLines] = size(omarLines);
    if numLines > 0
        fprintf("Fixing Path, Found %d lines starting with /home/omar\n", numLines);
        rmpath(omarLines{:});
    end

    addpath('/home/djhakim/gurobi811/linux64/matlab');
    addpath('/home/djhakim/lewisLab/geneProteinReaction/cobratoolbox');
    addpath('/home/djhakim/lewisLab/FBA');

    disp("Initializing Cobra Toolbox...");
    initCobraToolbox;

    disp("Reading sampleTPM.csv...");
    tpmData = readtable('sampleTPM.csv');
    cellLineNames = tpmData.Properties.VariableNames(3:end);
    [~,numCellLines] = size(cellLineNames);

    fprintf("Found %d cell lines.\n", numCellLines);

    disp("Loading iCHOv1.mat");
    load('iCHOv1.mat', 'iCHOv1');
    choModel = iCHOv1;

    disp("Loading ubiquity scores");
    ubiquityScore30 = readtable('ubiquityScore30.csv');

    disp("Constraining reaction bounds");
%     Overriding lb and ub for all reactions measured in Kallehauge study,
%     leaving all other lb and ub the same as the raw iCHOv1 model.
    choModel = setKallehaugeBounds(choModel);

    disp("Setting lower bound on biomass reaction");
    biomassRxn = find(choModel.c);
    choModel.lb(biomassRxn) = .0001; %#ok<FNDSB>    
    
    disp("Creating flux consistent model... (this usually takes a few minutes)");
    [fluxConsistentModel, fluxConsistentUbiquityScores] = createFluxConsistentModel(choModel, ubiquityScore30);

    disp("Creating tissue specific models");
    for i = cellLineStartIndex:cellLineEndIndex
        fluxConsistentUbiquityScore = fluxConsistentUbiquityScores{:,i};
        confidenceScore = zeros(size(fluxConsistentUbiquityScore));
        [tissueModel, coreRxn, nonCoreRxn,	zeroExpRxns, pruneTime, cRes] = mCADRE(fluxConsistentModel, fluxConsistentUbiquityScore, confidenceScore);

        results.tissueModel = tissueModel;
        results.coreRxn = coreRxn;
        results.nonCoreRxn = nonCoreRxn;
        results.zeroExpRxns = zeroExpRxns;
        results.pruneTime = pruneTime;
        results.cRes = cRes;

        filename = sprintf("model_%s.mat", cellLineNames{i});
        fprintf("Saving results to %s\n", filename);

        save(filename, "results");
    end
end
