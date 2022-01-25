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
// displs is displacement_vector
// rpp is rows_per_process
// np0 is number_of_process_with_rows=rpp
// np1 is number_of_process_with_rows=rpp+1
// var_perm is array to store temp_col_num
// rows_before is actual_posn in 1d view
// eff_row is proc_assign
// var_perm is perm_of_var

void perform_elimination(int recv_id, int id, int num_eq, float *proc_rows, float *proc_vals, int curr, float *recvd_row, int rows_per_proc, int num_proc, int *var_perm)
{

    swap(var_perm[(int)recvd_row[num_eq]], var_perm[recv_id]);
    ll i = 0;
    while (i < rows_per_proc)
    {
        if (num_proc * i != recv_id - id)
        {
            swap(proc_rows[(num_eq * i) + recv_id], proc_rows[(num_eq * i) + (int)recvd_row[num_eq]]);
            if (curr / num_proc > i)
            {
            }
            else
            {
                float piv_val = proc_rows[(num_eq * i) + recv_id];
                for (ll j = 0; j <= num_eq - 1; ++j)
                    proc_rows[j + (num_eq * i)] -= (recvd_row[j] * piv_val);
                proc_vals[i] = proc_vals[i] - recvd_row[num_eq + 1] * piv_val;
            }
        }

        i++;
    }
}

ll compute_pivot(int curr, int num_proc, int num_eq, float *proc_rows)
{
    float mx;
    ll row_id;
    row_id = curr / num_proc;
    mx = -1.0;
    ll pivot, i = curr;
    float val;
    while (i <= num_eq - 1)
    {
        val = proc_rows[i + (num_eq * row_id)];
        if (val > mx)
            pivot = i, mx = val;
        if (val <= 0)
            val = val * -1;
        i++;
    }
    return pivot;
}

void perform_division(int id, int curr, float *proc_rows, int pivot, int num_proc, int num_eq, int rows_per_proc, float *proc_vals)
{

    swap(proc_rows[pivot + ((curr / num_proc) * num_eq)], proc_rows[curr + ((curr / num_proc) * num_eq)]);
    ll i = curr;
    float piv_val = proc_rows[curr + ((curr / num_proc) * num_eq)];
    while (i <= num_eq - 1)
    {
        proc_rows[i + ((curr / num_proc) * num_eq)] /= piv_val;
        i++;
    }
    proc_vals[curr / num_proc] /= piv_val;
}

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
    // sending no. of equations to every process
    MPI_Bcast(&num_eq, 1, MPI_INT, root_rank, comm);

    // Defining variables for dividing processes
    ll rpp, rows_per_proc, np0, np1;
    ll var_perm[num_eq]; // stores the column no. in array
    for (ll i = 0; i < num_eq; ++i)
        var_perm[i] = i;
    rpp = num_eq / num_proc; // no of process that shoud be divided into processes
    // like 12/4=3

    // now basically what i m trying to do is pair conditions like
    // 15/4 so there will be like 3 processes with one more row and only
    // one process with 15/4

    ll rem_proc = num_eq % num_proc;
    np0 = num_proc - rem_proc;
    np1 = num_proc - np0;
    // if process rank is less than np1 so we should give it one more pocess
    rows_per_proc = ((id < np1) ? (1 + rpp) : rpp);
    ll divs[num_proc], displs[num_proc];
    memset(divs, 0, sizeof(divs));
    if (id = root_rank)
    {
        // now what i m trying to do is calc the final posn of row in linearized array
        memset(divs, 0, sizeof(divs));
        ll rows_before, eff_row;
        ll i = 0;
        displs[0] = 0; // used later for disp vector

        while (i < num_eq)
        {
            ll k = 0;
            rows_before = 0;
            if ((i % num_proc) > np1)
            {
                rows_before = np1 * (1 + rpp);
                rows_before = rpp * (i % num_proc) - rpp * np1 + rows_before;
            }
            else
                rows_before = (i % num_proc) * (1 + rpp);
            eff_row = divs[i % num_proc] + rows_before;
            divs[i % num_proc]++;
            // now we have assigned row in process , now we just read and allocate from imput
            while (k < num_eq)
            {
                fscanf(inputfp, "%f", &eq_mat[k + num_eq * eff_row]);
                k++;
            }
            fscanf(inputfp, "%f", &val_mat[eff_row]);
            i++;
            // now we have our A matrix and b matrix in desired form
        }
        // now we calculate the disp vector for scatterv
        ll itr = 1;
        while (itr < num_proc)
        {
            divs[itr - 1] = num_eq * divs[itr - 1];
            displs[itr] = divs[itr - 1] + displs[itr - 1];
            itr++;
        }
        divs[num_proc - 1] = num_eq * divs[num_proc - 1], fclose(inputfp);
    }

    // MPI_Barrier blocks all MPI processes in the given communicator until they all call this routine.
    // we are doing this so that our code input becomes in sync
    MPI_Barrier(comm);

    MPI_Finalize();
    return EXIT_SUCCESS;
}