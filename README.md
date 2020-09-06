# CGA-MWS-algorithm

# CGA-MWS algorithm JAVA source code

## A brief description of algorithm

* INPUT: A weighted non-binary mutation matrix or A mutation matrix, Parameter k.
* OUTPUT: Gene set of size k.

## Some features of the algorithm

* Integrate the two models (a new model based on the coefficient of variation and a original maximum weight submatrix solving model), and simply enter the set model name, which is globally applicable and easy to operate.
* Repeated operations of multiple groups are implemented in parallel through multiple threads to improve the efficiency of algorithm execution.
* Automatically execute p-value test after each algorithm run
* By adding the parameter of the number of times to execute the algorithm, the algorithm can be executed many times, and the results after each execution and the optimal value in the results of multiple execution algorithm can be printed.
* The genetic algorithm can be separated and executed independently, only in `GA_Algorithm.java` add the main function and input parameters.

## Format of TXT file for input algorithm
| Gene | TP53 | CDKN2A | CDKN2B| RB1 | CDK4| … |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| Sample_1 | 1 | 1.5 | 0 | 0 | 1 | … |
| Sample_2 | 0.45 | 1 | 0 | 0.4 | 1.5 | … |
| Sample_3 | 0 | 1.5 | 0 | 0.3 | 0 | … |
| … | … | … | … | … | … | … |
* If you want to use the original maximum weight submatrix model, enter a TXT file containing only `01` binary values similar to the table format above.

## Preparations before starting the program

1. The `main` method in `Run.java` is the entry to the whole program.

2. Enter the path to the TXT file at this location.
   ```Java
   String path = "A.txt;";
   ```
   
3. Parameter setting.

   Input parameters in this method.
   ```Java
   String[] paths = path.split(";");
	int g = 1126;
   int k = 2;
   int size = g / 2;
   r.run(paths, g, k, size, 500, 0.3, 10, 1000, "calfitness_Cov");   
   ``` 
   * The first   parameter:  The path of TXT file，
   * The second  parameter:  Number of genes in TXT file，
   * The third   parameter:  The size of Gene set (k)，
   * The fourth  parameter:  Population size (N, this value is the size of all populations combined)，
   * The fifth   parameter:  Iteration steps (maxg)，
   * The sixth   parameter:  Mutation probability (Pm)，
   * The seventh parameter:  Number of times the algorithm is executed，
   * The eighth  parameter:  Number of cycles when calculating p-value,
   * The ninth   parameter:  Model name ("calfitness_Cov", "calfitness_01"),
   * "calfitness_Cov": Model based on coefficient of variation,
   * "calfitness_01":  The original maximum weight sub-matrix solution model.
4. After setting the parameters, CGA-MWS algorithm can be executed.

## Some supplementary notes

* The project provides a TXT file `A.txt` containing only `01` values, which can be used as a test sample for algorithm testing.
* When testing is required, all `.java` files need to be downloaded.
* The CGA-MWS algorithm is based on JAVA8 implementation, please note the JAVA version when executing.
* The CGA-MWS algorithm printouts are shown below:

      NO.1 time 	The optimal gene set is：
      CDKN2B	CDK4	
      fitness:
      58.0	60.0	62.0	
      1 times, the average execution time is：0.21s
      p-value is: 1.0

