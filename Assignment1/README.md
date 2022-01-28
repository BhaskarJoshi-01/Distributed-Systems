## DS ASSIGNMENT 1
### Bhaskar Joshi
### 2019111002

#### Brief report

#### Instruction to run the code
Firstly install OpenMPI 
The following link can be helpful 
https://edu.itp.phys.ethz.ch/hs12/programming_techniques/openmpi.pdf

Considering MPI is installed
- First Question can be compiled using:

     ```mpic++ 2019111002_1.cpp && mpirun -np "PLACE DESIRED No. OF PROCESSES" --oversubscribe ./a.out "INPUT FILE PATH" "OUTPUT FILE PATH"```

- Second Question can be compiled using:

     ```mpic++ 2019111002_2.cpp && mpirun -np "PLACE DESIRED No. OF  PROCESSES" --oversubscribe ./a.out "INPUT FILE PATH" "OUTPUT FILE PATH"```

- Similarly Third Question can be compiled using:

     ```mpic++ 2019111002_3.cpp && mpirun -np "PLACE DESIRED No. OF  PROCESSES" --oversubscribe ./a.out "INPUT FILE PATH" "OUTPUT FILE PATH"```

The above codes can be executed considering the user has entered the correct arguments and correct path.

#### Question 1 Brief 
We were supposed to write a program to check if a number is prime. <br>
Let the given number be N.<br>
Then I have used `MPI_Bcast()` to brodcast value of N to all processes and then iterated through sqrt(N) and then divided the number possible ranges [start, end] to process according to its rank and check if there is any number in the range that divides sqrt(N) if there is the flag is raised to 1.<br>
Now I have used `MPI_Reduce()` so that if reduction result is non zero the number has to be non prime else the number will be prime.

#### Question 2 Brief
In this question we were given an undirected binary weighted graph G, find the number of cliques (of all possible weights) in the graph with sizes 3 and 4.<br>
After getting n and e as input, I have used `MPI_Bcast()` to broadcast value of n to all processes. Next thing was to generate the Adjacency Matrix and broadcast it to all the processes. Now every process has n and adj matrix.<br>
Now I have done the distributed work on the basis of rank to all the processes, and then obtained cliques for vertices using and it has all the necessary data required for execution.<br>

```get_cliques_for_vertices(start_vertex, end_vertex, n, adjacency_matrix, clique_count);```

``` void combinations(int arr[], int data[], int start, int end, int index, int r, int adjacency_matrix[], int clique_count[], int n)```<br>
This function calls combination function which returns the possible number of cilques of size 3 and 4 via backtracking.<br>
Then the result is gathered using `MPI_Reduce()` and outputed in the output file.

#### Question 3 Brief 
In this question we were expected to write a program to find the solution to the set of linear equations Ax = b using the  Gaussian Elimination algorithm followed by back substitution.<br>
First we take n from the input file and broadcast this using `MPI_Bcast()` to all the processes then we linearize other input i.e A matriz and B matrix but prior to this we first assign P_0 and P_1 it is the number of process with n/np size and n/np+1 respectevely and then we assign our rows based on rank and also store the relative position of the column oof current row. <br>  
Then for root rank zero we read the inputs and store them into flattend array defined above, the advantage of this is that now our processes are in cyclic order rather than being chunks of contiguous blocks and this could prevent one process to work after getting immediately previous input rather than waiting for the entire chunk of rows to get the task done. Also we define the displacement array the use of this will be in `MPI_Scatterv()`.Basically what it does is that it sends three different kind of values in array combined i.e value of row (n entities) , pivot (1 entity), corresponding b value (1 entity) in the exact same order. We then scatterv the value of A and b matrices according to size and displacement array already defined above. <br>
The next part is to generation of upper triangular matrix now after the communication described , the data is communicated to immediate next process by rank.<br>
as a row is recieved to a process it first communicates it to the next row and then preforms the elimination whose indices are greater than the indices of recieved rows. Once it is determined that row of current process can be used for elimination of other rows , it places pivot then divides entire row by the pivot. This information is then communicated forward. `MPI_Send` and `MPI_Recv` are used in the process. <br>
Finally we have back substitution and it is smilar to the previous part except the fact that it is in reverse order only pivot values are transmitted upward instead of rows and the values of the variables are determined using this procedure hence the system is also consistent and satisfied.
For this apart from `MPI_Send` and `MPI_Recv` we also use `MPI_Gatherv` for getting the results.
