function [model, ubiquityScore, confidenceScores] = runStanDep(numClusters)
    
    [expData, model] = loadForStanDep();
    
    %Flatten impossible numbers, sometimes we see TPM data of 1xe-300.
    %This obviously didn't come from a real measurement.  Removing this
    %type of data point prevents issues later in the workflow when
    %calculating standard deviations of gene expression in the logarithmic
    %frame.  
    
    badTPM = expData.valuebyTissue > 0 & expData.valuebyTissue < 0.0000001;
    expData.valuebyTissue(badTPM) = 0;

    modelData = getModelData(expData,model);
    
    spec = getSpecialistEnzymes(model);  
    prom = getPromEnzymes(model);
    enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

    edgeX = [-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6]; % bins  
    distMethod = 'euclidean'; % distance method  
    linkageMethod = 'complete'; % linkage metric for hierarchical clustering
    
    clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,numClusters,distMethod,linkageMethod);

    
    coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,"./output",[],true,0,[1 1]); 
    calcJaccardSimilarity(coreRxnMat,enzymeData.Tissue,'matrix',true);
    ubiquityScore = getUbiquityScore(clustObj,edgeX,model); % calculate ubiquity score
    confidenceScores = zeros(size(ubiquityScore));
    
    %Update the ubiquity score of the biomass reaction(s) to be 1, standep
    %can't recognize this reaction because it has no genes associated with
    %it.    
    ubiquityScore(model.c == 1, :) = 1;
    
    %Force enable reactions r1664 to r2071 -- These are all dependent on
    %the same two genes, but by default one of these genes is split by
    %StanDep in a way we disagree with when looking at the distributions
    %manually.  
    forceEnable = model.rxns >= "r1664" & model.rxns <= "r2071";
    ubiquityScore(forceEnable, :) = 1;
    
    %TODO FIXME HACK: we only noticed the r1664 to r2071 issue because it
    %dramatically skewed our output results.  We need to somehow validate
    %that this type of issue isn't prevalent in our results.  
end