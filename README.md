# stochastic_nonlinear_chaos
Implementation of some nonlinear equation systems using stochastic computing, in Matlab

For more information, have a look at our paper **Stochastic computing implementation of chaotic systems** (https://doi.org/10.3390/math9040375). If you like this code and use it, please cite the original source: 
    Camps, O.; Stavrinides, S.G.; Picos, R. Stochastic Computing Implementation of Chaotic Systems. Mathematics 2021, 9, 375. 


- The lorenz/ directory contains the implementation of the Lorenz system. There are two files there: the lorenz.m and the lorenz_s.m. The first one implements the lorenz system "as it is", with the option of setting a specific number of bits for all the operations. The lorenz_s.m implements the same system using stochastic computing. Notice that all the operations are implemented, not as operators, but as functions, so the actual integration is done over the same formal functions. The stochastic implementations also write the output to a file and, if interrupted and restarted, they try to read this file and continue from the last saved point. Notice that, since we're forcing a given precision to the numbers, no actual information is lost.

- The shi-mo/ directory contains the implementation of the Shimizu-Morioka system (Shimizu, T.; Morioka, N. On the bifurcation of a symmetric limit cycle to an asymmetric one in a simple model. Phys. Lett. A 1980, 76, 201–204. [Google Scholar] [CrossRef]) , using the same conventions as for the Lorenz, including the auto-saving and recovering feature.

- The zhang-tang/ directory contains the implementation of the Zhang-Tang system (scitation pending), using the same conventions and implementation than for the Lorenz, including the auto-saving and recovering feature.
