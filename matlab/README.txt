In order to run the matlab code you must:

	Retrieve the ICHOv1.mat file representing the CHO genome scale metabolic model from the BiGG database

	Install the cobra toolbox and add it to your matlab path
		Install any of the necessary cobra toolbox mathematical solvers, recommend gurobi
	
	Install StanDep and add it to your matlab path

	Acquire RNA Seq data (Table4_omics.xlsx)

You can then
	runStanDep (modifying number of clusters and bin edges as needed) to convert RNA seq data into ubiquity scores used by mCADRE.  
	generateCellLineModels to use ubiquity scores from StanDep and reaction bounds from Kallehauge to generate context specific models with mCADRE.  

Once you have these models, you can print them to csv files with 
	printModels

And analyze them with the code in the python folder.