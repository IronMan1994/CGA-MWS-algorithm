# CGA-MWS-algorithm

# CGA-MWS algorithm JAVA source code

## A brief description of algorithm

* INPUT: A weighted non-binary mutation matrix or A mutation matrix.
* OUTPUT: Gene set of size k.

## Some features of the algorithm

* Integrate the two models (a new model based on the coefficient of variation and a original maximum weight submatrix solving model), and simply enter the set model name, which is globally applicable and easy to operate.
* Repeated operations of multiple groups are implemented in parallel through multiple threads to improve the efficiency of algorithm execution.
* The genetic algorithm can be separated and executed independently, only in `GA_Algorithm.java` add the main function and input parameters.
* By adding the parameter of the number of times to execute the algorithm, the algorithm can be executed many times, and the results after each execution and the optimal value in the results of multiple execution algorithm can be printed.

## Format of TXT file for input algorithm
| Gene | TP53 | CDKN2A | CDKN2B| RB1 | CDK4| … |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| Sample_1 | 1 | 1.5 | 0 | 0 | 1 | … |
| Sample_2 | 0.45 | 1 | 0 | 0.4 | 1.5 | … |
| Sample_3 | 0 | 1.5 | 0 | 0.3 | 0 | … |
| … | … | … | … | … | … | … |


'You can use a TXT file in the same format as the sample file provided for code testing.
Or if you just want to run the code for testing, you can use the provided TXT file directly.
The "main" function in Run.java is the entry point of the entire program. 
Please adjust the parameters in the "main" function.'
