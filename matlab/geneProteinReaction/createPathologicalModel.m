function [badModel] = createPathologicalModel()
    badModel.mets = cell(4,1);
    badModel.mets{1,1} = 'A';
    badModel.mets{2,1} = 'B';
    badModel.mets{3,1} = 'C';
    badModel.mets{4,1} = 'D';
    
    badModel.metNames = cell(4,1);
    badModel.metNames{1,1} = 'A';
    badModel.metNames{2,1} = 'B';
    badModel.metNames{3,1} = 'C';
    badModel.metNames{4,1} = 'D-Biomass';
    
    badModel.metFormulas = cell(4,1);
    badModel.metFormulas{1,1} = 'C2';
    badModel.metFormulas{2,1} = 'C2';
    badModel.metFormulas{3,1} = 'C2';
    badModel.metFormulas{4,1} = 'C2';

    badModel.metCharge = cell(4,1);
    badModel.metCharge{1,1} = 0;
    badModel.metCharge{2,1} = 0;
    badModel.metCharge{3,1} = 0;
    badModel.metCharge{4,1} = 0;
    
    badModel.genes = cell(4,1);
    badModel.genes{1,1} = 'g1';
    badModel.genes{2,1} = 'g2';
    badModel.genes{3,1} = 'g3';
    badModel.genes{4,1} = 'g4';
    
    badModel.rxnGeneMat = zeros(6,4);
    badModel.rxnGeneMat(1,1) = 1;
    badModel.rxnGeneMat(2,2) = 1;
    badModel.rxnGeneMat(3,3) = 1;
    badModel.rxnGeneMat(6,4) = 1;
    
    badModel.rules = cell(6, 1);
    badModel.rules{1,1} = 'x(1)';
    badModel.rules{2,1} = 'x(2)';
    badModel.rules{3,1} = 'x(3)';
    badModel.rules{4,1} = '';
    badModel.rules{5,1} = '';
    badModel.rules{6,1} = 'x(4)';

    badModel.grRules = cell(6,1);
    badModel.grRules{1,1} = 'g1';
    badModel.grRules{2,1} = 'g2';
    badModel.grRules{3,1} = 'g3';
    badModel.grRules{4,1} = '';
    badModel.grRules{5,1} = '';
    badModel.grRules{6,1} = 'g4';
    
    badModel.rxns = cell(6,1);
    badModel.rxns{1,1} = 'rAIn';
    badModel.rxns{2,1} = 'rAToB';
    badModel.rxns{3,1} = 'rAToC';
    badModel.rxns{4,1} = 'rBOut';
    badModel.rxns{5,1} = 'rCOut';
    badModel.rxns{6,1} = 'rBCToD';
    
    badModel.rxnNames = cell(6,1);
    badModel.rxnNames{1,1} = 'rnAIn';
    badModel.rxnNames{2,1} = 'rnAToB';
    badModel.rxnNames{3,1} = 'rnAToC';
    badModel.rxnNames{4,1} = 'rnBOut';
    badModel.rxnNames{5,1} = 'rnCOut';
    badModel.rxnNames{6,1} = 'rnBCToD';
    
    badModel.subSystems = cell(6,1);
    badModel.subSystems{1,1} = 'DAN';
    badModel.subSystems{2,1} = 'DAN';
    badModel.subSystems{3,1} = 'DAN';
    badModel.subSystems{4,1} = 'DAN';
    badModel.subSystems{5,1} = 'DAN';
    badModel.subSystems{6,1} = 'DAN';
    
    badModel.csense = transpose('EEEE');
    
    badModel.S = zeros(4,6);
    badModel.S(:,1) = [1,0,0,0];
    badModel.S(:,2) = [-1,1,0,0];
    badModel.S(:,3) = [-1,0,1,0];
    badModel.S(:,4) = [0,-1,0,0];
    badModel.S(:,5) = [0,0,-1,0];
    badModel.S(:,6) = [0,-1,-1,0];
    
    badModel.lb = zeros(6,1);
    badModel.lb(1) = 100;
    badModel.lb(2) = 0;
    badModel.lb(3) = 0;
    badModel.lb(4) = 0;
    badModel.lb(5) = 0;
    badModel.lb(6) = 0;
    
    badModel.ub = zeros(6,1);
    badModel.ub(1) = 100;
    badModel.ub(2) = 1000;
    badModel.ub(3) = 1000;
    badModel.ub(4) = 1000;
    badModel.ub(5) = 1000;
    badModel.ub(6) = 1000;
    
    badModel.b = zeros(4,1);
    
    badModel.c = zeros(6,1);
    badModel.c(6) = 1;
    
    badModel.rev = zeros(6,1);
    badModel.description = 'PathologicalmCADRETest';
    
    results = zeros(100,1);
    
    for i = 1:100    
        ubRandom = zeros(6,1);
        ubRandom(1,1) = rand;
        ubRandom(2,1) = rand;
        ubRandom(3,1) = rand;
        ubRandom(4,1) = -1;
        ubRandom(5,1) = -1;
        ubRandom(6,1) = rand;
        
        confidence = zeros(6,1);
        
        contextSpecificModel = mCADRE(badModel, ubRandom, confidence);
        fba = optimizeCbModel(contextSpecificModel);
        results(i,1) = fba.f;
    end
    for i = 1:100
        fprintf("Optimal Biomass: " + results(i,1) + "\n");
    end
    fprintf("Total mCADRE Biomass: " + sum(results) + "\n");
    fba = optimizeCbModel(badModel);
    fprintf("RAW Optimal Biomass: " + fba.f + "\n");
end 