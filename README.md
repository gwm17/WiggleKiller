# WiggleKiller
WiggleKiller is a program designed to address differential non-linearity in the SESPS focal plane detector position spectrum (henceforth and forevermore known as the dreaded wiggles). The priciple is that the DNL can be found as
deviations in the PSD vs. time for each delay line signal. This can then be calibrated away using the following relation:

t' = t + lambda*psd

This essentially creates a re-calibrated time for each signal. Lambda is a free parameter to be deterimed either through qualitative testing or using the given optimization method (with a warning that the optimization might still
have quirks). Once a good lambda is deterimed, in can in principle be applied to several data sets using 

t' = t + dt

where dt is the constant lambda*psd. This is not in general applicable, but for data sets with similar particles at similar energies (read: similar pulse shapes), it should be fine to use. However one could also just use the parameter
optimization/determination to find a new calibration. The final product can either be a file of histograms (for tests) or a new cleaned-up ROOT tree in a file.

Data is assumed to be of the type given from the GWM_EventBuilder SPS program. Additionally, cuts will be applied to the data, and the final calibrated data file will have these cuts applied (i.e. permanently cut on particle ID).  

## Optimization method
Optimization is done using a simulated annealing method. The optimization parameter is the peak-to-trough height of one of the wiggles (but better parameters ar being tested). The peak and trough location are manually entered into the
code in the WiggleKiller.cpp file. Simulated annealing is slower than other methods, however it requires very little knowledge/assumptions about the shape of the function to be minimized, and does a better job approaching a global
minimum than other simple methods such as a Golden-Section Search.


## Building and Execution
WiggleKiller uses CMake. To build WiggleKiller run the following from the repository directory:

```
mkdir build
cd build
cmake ..
make
```
To run WiggleKiller use the following command structure:

./bin/WiggleKiller `--<option>` input.txt

`--<option>` is one of the following options:

1. kill: runs the single-shot correction of the lambda values in the input file. 
2. optimize: runs the optimization using the initial guess of the lambda values in the input file. 
3. clean: generates a cleaned-up data file with a ROOT tree using the calibration files specified in the input file.
4. clean-no-tree: generated a cleaned-up data file with only histograms (no ROOT tree) using the calibration files specified in the input file. 

An example input file is included. The input specifies directory locations for input data and output data.
