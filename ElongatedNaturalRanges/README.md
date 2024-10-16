# ConcavityControlsComplexity
# In ElongatedNaruralRanges 

1. A Matlab format Supplementary Table 1, containing the data on the 18 elongated mountain ranges, is stored in files:
LinearRanges145.mat and LinearRanges145WithAI.mat. These files can be used to reproduce figures 2a and 6, using the matlab script: PlotNaturalRangesFigs.


2. The matlab function ConcavityAndComplexityBasedOnPartial.m generates the data that was used to feed the rows (each row per mountain range) in the table above. 
An input to the function is a mat file whose name has the format: RangenameDataPartial145.mat, where Rangename is the name of the range. 
The function CallConcavityComplexityPartial.m contains code lines that run ConcavityAndComplexityBasedOnPartial.

3. The input files RangenameDataPartial145.mat for the function ConcavityAndComplexityBasedOnPartial.m are included in the current branch. Each such file corresponds to a single range. These files were generated by interacting with the relevant DEM and choosing the main basins for the analysis. The input files corresponding to the (i) Central Mountain Range, Taiwan, (ii) Finisterre Range, Papua New Guinea, and (iii) Sierra Madre del Sur, Mexico are too large to be hosted in GitHub and can be requested from the corresponding author.

4. Other matlab functions included in the current branch are service functions used by ConcavityAndComplexityBasedOnPartial.m

5. The function PlotNiceRanges.m plots the drainage network analyzed in each range over the range relief map.


The codes run on Topotoolbox environment. 
