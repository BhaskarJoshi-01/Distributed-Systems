#include <mpi.h>
#include <bits/stdc++.h>
#include <fstream>
#define ll long long int

// num_proc is the number_of_processes
// id is the rank of process
// inputfp is the input_file_pointer
// eq_mat is the A_matrix
// val_mat is the b_matrix
// divs is no_of_coefficients_per_process
// displs is begining_index_of_coeff_of_each_process
// rpp is rows_per_process
// np0 is number_of_process_with_rows=rpp
// np0 is number_of_process_with_rows=rpp+1

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    int num_proc, id;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &num_proc);
    MPI_Comm_rank(comm, &id);
    int root_rank = 0;
    // Reading input from file
    FILE *inputfp = NULL;
    ll num_eq;
    // reading number of equations from file first and creating
    // variables accordingly for taking matrix as input and distributing

    if (id == root_rank)
    {
        inputfp = fopen(argv[1], "r");
        fscanf(inputfp, "%lld", &num_eq);
    }

    ll oned_size = pow(num_eq, 2);
    // converting n*n A matrix that was input to 1*n matrix

    float eq_mat[oned_size], val_mat[num_eq];
    ll divs[num_proc], displs[num_proc];
    // sending no. of equations to every process
    MPI_Bcast(&num_eq, 1, MPI_INT, root_rank, comm);
    
    //Defining variables for dividing processes
    ll rpp,np0,np1;
    rpp = num_eq / num_proc; // no of process that shoud be divided into processes
    // like 12/4=3 

    np0 = num_proc - (num_eq % num_proc); // Number of processes with rows = rpp (np0),
    np1 = num_proc - np0;


    MPI_Finalize();
    return EXIT_SUCCESS;
}