# ConcavityControlsComplexity
# In GreedyOCN


The greedy version of Optimal Channel Network Generator with the implementation of delta chi and delta L caculation. 
For a quick and easy run, try: 

[E_all,Prob_all,drca_all,theta] = runSimulatedAnnealing(80,30,1,300,1e9,0.9);

To generate data similar to that used in the manuscript try: 
[E_all,Prob_all,drca_all,theta] = runSimulatedAnnealing(200,60,1,300,1e9,0.9);

The data used in the manuscript is stored in files:
SA_160423_010203.mat
SA_160423_040506.mat
SA_160423_070809.mat

to produce the figures from the manuscript, run PlotResults.m

The codes run on Topotoolbox environment 
