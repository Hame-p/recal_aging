# recal_aging
This project contains codes and data from the preprint:
Sensory- and memory-related drivers for altered ventriloquism effects and aftereffects in older adults
Hame Park, Julia Nannt, Christoph Kayser
doi: https://doi.org/10.1101/2020.02.12.945949

It contains two folders, Code and Data.
Code contains the codes to do the GLMM and plot figures 2-4 in the paper.
Data contains behavioral data and the result data (parameters, negative log-likelihoods) of the causal inference model
and its simulation. 

Data:
1. aging.mat
This is the behavioral data from the young and older adults performing a rapid ventriloquism effect/aftereffect task.
All stimulus location and responses, biases VE and VAE are in visual degrees.
0 is center, negative values are to the left, positive values are to the right.
For basic paradigm, see https://elifesciences.org/articles/47001, Figure 1 & Methods

	1-1. The dataset is a cell array of two sub-datasets.
	The lowest level of array is the subject data.
	N sub-dataset = 2
	sub-dataset1: N subject = 22, healthy young adults
	sub-dataset2: N subject = 21, healthy older adults

	1-2. Each subject data is a matrix and the first 8 columns are:
	c1: AV trial V stimulus location
	c2: AV trial A stimulus location
	c3: A trial A stimulus location
	c4: dVA: AV V - AV A (a.k.a, audio-visual discrepancy)
	c5: AV trial response
	c6: VE
	c7: A trial response
	c8: VAE

2. VE_BCI_result.mat
This file contains the model fitting result.

	2-1. AICc: corrected Akaike information criterion 
	(2 age groups x 3 models (1.model averaging, 2.model selection, 3.probability matching))
	
	2-2. BIC: Bayesian information criterion (same format as AICc)
	
	2-3. bvp: best model parameters 
	(2 age groups x 5 parameters: Sigma_A, Sigma_V, Sigma_P, Mu_P, Pcom, best model index (1-3 from 2-1))
	
	2-4. mNegLL: mean negative log-likehoods (same format as AICc)
	
	2-5. mParams: mean model parameters 
	(2 age groups x 3 models x 5 parameters: Sigma_A, Sigma_V, Sigma_P, Mu_P, Pcom, best model index (1-3 from 2-1))
	
	2-6. out, post: output from VBA_groupBMC.m (VBA-toolbox)

3. VE_simSD.mat
This file contains simulated VE with fitted parameters from the BCI model.

	3-1. xpsd: simulated VE for 9 audio-visual discrepancies. 
	2 age groups x 3 models, contains matrix of 9 discrepancies x (mean , SD) x subjects
	
	3-2. bxpsd: same as xpsd but with the best fitted parameters for each subject. 

4. nobs_A_AVtrial.mat
nobs: number of observations for each subject used to fit the model. 
