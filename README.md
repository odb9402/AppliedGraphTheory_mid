# Finding the Shortest Cycle in Delaunay Graph with Constraint

A mid-term project of PNU graduated course "Applied graph theory" by prof. Cho



## Problem

Let assume that there is a set of two-dimensional points. Delaunay Graph(DG) is defined by using those geometrical points. If there is a cycle on DG, *dist* is defined as the largest value of every shortest path between the cycle and a vertex of DG. 

We have to minimize the length of a cycle on DG with a constraint that *dist* should be shorter than a constraint value. In this project, we consider this problem as Integer Programming with Constraint. We applied Genetic Algorithm to solve the optimization problem.



## Method





## Short description to test

Input file and output file follow the format of the course.

#### Python3 dependency
- deap
- networkx
- numpy
- scipy
- matplotlib

`pip install deap networkx numpy scipy matplotlib`

#### Set input file name from the source

`P, threshold = read_points("Input_file_name")`



#### Check the output

`write_cycle(cycle, "output_file_name")`



#### Adjust generations for the Genetic Algorithm

In this source, there are several parameters to adjust. In the report for the submission, all parameter was fixed as below ( same with the source ). During the test for the algorithm, pop_size and gen_num can be adjusted.



`pin_num = 30`

`pop_size = 300`

`gen_num = 500`

