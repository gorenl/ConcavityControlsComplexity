# ConcavityControlsComplexity

1. DAC simulations data is stored in file Dac_Data_Apr_2023.mat. To plot run DelLDelchi.

2. To regenerate the data, use the function callDeltaChiDeltaLToShare that uses the folders Theta0.*ASCII. These folders contain  simulations output ASCII files. 

3. To rerun the simulations, we supply the source code. To change the value of theta, change line 27 in file initialize_parameters.f90 to the desired m value. Check line 49 for the required erodibility (see supplementary table S2 to know the value), and save the file. Compile the code (run makes (requires gfortran) and run it with (./DivideAndCapture.exe). Simulation output stored in folder RUN5 can be viewed in Paraview  
