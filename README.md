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
All stimulus location and responses, biases VE and VAE are in visual degrees.
0 is center, negative values are to the left, positive values are to the right.
For basic paradigm, see https://elifesciences.org/articles/47001, Figure 1 & Methods

1. The dataset is a cell array of two sub-datasets.
The lowest level of array is the subject data.
N sub-dataset = 2
sub-dataset1: N subject = 22, healthy young adults
sub-dataset2: N subject = 21, healthy older adults

2. Each subject data is a matrix and the first 8 columns are:
c1: AV trial V stimulus location
c2: AV trial A stimulus location
c3: A trial A stimulus location
c4: dVA: AV V - AV A (a.k.a, audio-visual discrepancy)
c5: AV trial response
c6: VE
c7: A trial response
c8: VAE

