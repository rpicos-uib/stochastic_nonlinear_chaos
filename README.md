# stochastic_nonlinear_chaos
Implementation of some nonlinear equation systems using stochastic computing, in Matlab

For more information, have a look at our paper **Stochastic computing implementation of chaotic systems** (https://doi.org/10.3390/math9040375). If you like this code and use it, please cite the original source: 
    Camps, O.; Stavrinides, S.G.; Picos, R. Stochastic Computing Implementation of Chaotic Systems. Mathematics 2021, 9, 375. 


- The lorenz/ directory contains the implementation of the Lorenz system. There are two files there: the lorenz.m and the lorenz_s.m. The first one implements the lorenz system "as it is", with the option of setting a specific number of bits for all the operations. The lorenz_s.m implements the same system using stochastic computing. Notice that all the operations are implemented, not as operators, but as functions, so the actual integration is done over the same formal functions.

- The shi-mo/ directory contains the implementation of the Shimizu-Morioka system, using the same conventions as for the Lorenz


