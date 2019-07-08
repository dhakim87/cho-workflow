function [tissueModel, coreRxn, nonCoreRxn,	zeroExpRxns, pruneTime, cRes] = runMCADRE(fluxConsistentModel, fluxConsistentUbiquityScores, index)

%Just run mCADRE for the first cell line for now.
ubiquityScore = fluxConsistentUbiquityScores{:,index};
confidenceScore = zeros(size(ubiquityScore));
[tissueModel, coreRxn, nonCoreRxn,	zeroExpRxns, pruneTime, cRes] = mCADRE(fluxConsistentModel, ubiquityScore, confidenceScore);

end