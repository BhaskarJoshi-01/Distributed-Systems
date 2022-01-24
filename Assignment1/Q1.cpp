#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
#define ld long double
// This will work for every those number in the range of int
// This function is to check if a number is prime or not
// in range start->end , if it is prime it returns 1 else 0

int isPrime(int start, int end, int num)
{
    int flag = 0; //prime
    if (num == 1)
    {
        return 1; // non prime : corner case
    }
    for (int i = start; i <= end; i++)
    {
        if (num % i == 0 && num != 1 && i != 1)
        {
            flag = 1; // non prime
            // cout << start << " " << end << num << " " << i << endl;
            break;
        }
    }
    return flag;
}

int main(int argc, char *argv[])
{
    // Initialize MPI

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int numprocs, rank;
    // get size of the current communicator
    MPI_Comm_size(comm, &numprocs);
    // get current process rank
    MPI_Comm_rank(comm, &rank);

    // handeling corner cases
    if ( numprocs >= 12)
    {
        printf("Incorrect Input format \nGive only input as: \n mpirun -np 11 ./a.out <input file> <output file>. \n");
        MPI_Abort(comm, EXIT_FAILURE);
    }
    // Now taking inputs in rootrank i.e 0 and broadcasting it
    // need to edit these informations
    int N;
    if (rank == 0)
    {
        // cin >> N;
        fstream file;
        string line;
        file.open(argv[1]);
        if (file.is_open())
        {
            getline(file, line);
        }
        N = stoi(line);
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //The highest value of N can be 1e9 as given in question
    //so for general purpose I have not taken sqrt(1e9) i.e 1e3
    int sqrt_N = sqrt((ld)N);
    int a = rank * (sqrt_N / numprocs) + 1;
    if (rank == 0)
    {
        a = 2;
    }
    int b = (rank + 1) * (sqrt_N / numprocs);

    if (rank == numprocs - 1)
    {
        b = sqrt_N;
    }
    int my_result = isPrime(a, b, N);
    int reduction_result;
    // cout << my_result << " " << rank << endl;
    MPI_Reduce(&my_result, &reduction_result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {

        freopen(argv[2], "w", stdout);
        cout << ((reduction_result > 0) ? "NO" : "YES") << endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS; // defined as 0 in stdlib
}