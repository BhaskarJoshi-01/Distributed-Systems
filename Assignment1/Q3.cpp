#include <mpi.h>
#include <bits/stdc++.h>
#include <fstream>
#define ll long long int
#define ld double
using namespace std;

auto Eliminate = [](int recv_rank, int rank, int no_of_eq, ld *proc_rows, ld *proc_vals, int curr, ld *recvd_row, int rows_per_proc, int no_of_processes, int *permut_variables)
{
    swap(permut_variables[(int)recvd_row[no_of_eq]], permut_variables[recv_rank]);
    ll i = 0;
    while (i < rows_per_proc)
    {
        if (no_of_processes * i != recv_rank - rank)
        {
            swap(proc_rows[(no_of_eq * i) + recv_rank], proc_rows[(no_of_eq * i) + (int)recvd_row[no_of_eq]]);
            if (curr / no_of_processes > i)
            {
            }
            else
            {
                ld piv_val = proc_rows[(no_of_eq * i) + recv_rank];
                for (ll j = 0; j <= no_of_eq - 1; ++j)
                    proc_rows[j + (no_of_eq * i)] -= (recvd_row[j] * piv_val);
                proc_vals[i] = proc_vals[i] - recvd_row[no_of_eq + 1] * piv_val;
            }
        }
        i++;
    }
};

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    int no_of_processes, rank;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &no_of_processes);
    MPI_Comm_rank(comm, &rank);
    int root_rank = 0;
    // Reading input from file
    FILE *read_inp = NULL;
    ll no_of_eq = 0;
    // reading number of equations from file first and creating
    // variables accordingly for taking matrix as input and distributing

    if (rank == root_rank)
    {
        read_inp = fopen(argv[1], "r");
        fscanf(read_inp, "%lld", &no_of_eq);
    }
    // printf("First4");

    ll oned_size = pow(no_of_eq, 2);
    // converting n*n A matrix that was input to 1*n^2 matrix

    ld A_matrix[oned_size], b_matrix[no_of_eq];
    // sending no. of equations to every process
    MPI_Bcast(&no_of_eq, 1, MPI_INT, root_rank, comm);
    // printf("First5 %d\n", rank);

    // Defining variables for dividing processes
    ll rows_per_process, rows_per_proc, P_0, P_1;
    int permut_variables[no_of_eq]; // stores the column no. in array
    for (ll i = 0; i < no_of_eq; ++i)
        permut_variables[i] = i;
    rows_per_process = no_of_eq / no_of_processes; // no of process that shoud be divided into processes
    // like 12/4=3

    // now basically what i m trying to do is pair conditions like
    // 15/4 so there will be like 3 processes with one more row and only
    // one process with 15/4

    ll A_SZ, B_SZ, REC_SZ, rem_proc = no_of_eq % no_of_processes;
    P_0 = no_of_processes - rem_proc;
    P_1 = no_of_processes - P_0;
    // if process rank is less than P_1 so we should give it one more pocess
    rows_per_proc = ((rank < P_1) ? (1 + rows_per_process) : rows_per_process);
    // fprintf(stderr, "hello %d %lld\n", rank, rows_per_proc);
    int chunks_per_proc[no_of_processes], displacement_vector[no_of_processes];
    memset(chunks_per_proc, 0, sizeof(chunks_per_proc));
    if (rank == root_rank)
    {
        // now what i m trying to do is calc the final posn of row in linearized array
        memset(chunks_per_proc, 0, sizeof(chunks_per_proc));
        ll actual_posn_1d, proc_assign;
        ll i = 0;
        displacement_vector[0] = 0; // used later for disp vector
        while (i < no_of_eq)
        {
            ll k = 0;
            actual_posn_1d = 0;
            if ((i % no_of_processes) > P_1)
            {
                actual_posn_1d = P_1 * (1 + rows_per_process);
                actual_posn_1d = rows_per_process * (i % no_of_processes) - rows_per_process * P_1 + actual_posn_1d;
            }
            else
                actual_posn_1d = (i % no_of_processes) * (1 + rows_per_process);
            proc_assign = chunks_per_proc[i % no_of_processes] + actual_posn_1d;
            chunks_per_proc[i % no_of_processes]++;
            // cerr << "\ndivs update" << chunks_per_proc[i % no_of_processes] << endl;
            // now we have assigned row in process , now we just read and allocate from imput
            while (k < no_of_eq)
            {
                fscanf(read_inp, "%lf", &A_matrix[k + no_of_eq * proc_assign]);
                k++;
            }
            fscanf(read_inp, "%lf", &b_matrix[proc_assign]);
            i++;
            // now we have our A matrix and b matrix in desired form
        }
        // now we calculate the disp vector for scatterv
        ll itr = 1;
        while (itr < no_of_processes)
        {
            chunks_per_proc[itr - 1] = no_of_eq * chunks_per_proc[itr - 1];
            displacement_vector[itr] = chunks_per_proc[itr - 1] + displacement_vector[itr - 1];
            itr++;
        }
        chunks_per_proc[no_of_processes - 1] = no_of_eq * chunks_per_proc[no_of_processes - 1], fclose(read_inp);
    }

    // MPI_Barrier blocks all MPI processes in the given communicator until they all call this routine.
    // we are doing this so that our code input becomes in sync
    MPI_Barrier(comm);
    // next we have distributed the A matrix to other porcesses
    A_SZ = rows_per_proc * no_of_eq;
    B_SZ = A_SZ / no_of_eq;
    ld proc_rows[A_SZ]; // part of A matrix for a process
    ld proc_vals[B_SZ]; // part of B matrix for a process

    // if (rank == 0)
    // {
    //     // for (int i = 0; i < A_SZ; i++)
    //     // {
    //     //     cerr << A_matrix[i] << " ";
    //     // }
    //     // cerr << "\n";
    //     // cerr << "\n"
    //         //  << rank << "DIVS: ";
    //     for (int i = 0; i < no_of_processes; i++)
    //     {
    //         cerr << chunks_per_proc[i] << " ";
    //     }
    //     cerr << endl;
    //     cerr << "disps\n";
    //     for (int i = 0; i < no_of_processes; i++)
    //     {
    //         cerr << displacement_vector[i] << " ";
    //     }
    //     cerr << "\n";
    // }

    // MPI_Scatterv is in which the data dispatched from the root process
    // can vary in the number of elements, and the location from which load
    // these elements in the root process buffer.
    MPI_Scatterv(A_matrix, chunks_per_proc, displacement_vector, MPI_DOUBLE, proc_rows, A_SZ, MPI_DOUBLE, root_rank, comm);
    // cerr << "hi " << rank << endl;

    // similarly we need to do for B matrix
    // as chunks_per_proc and disps were multiplied by no_of_eq above
    // changing them back
    for (ll i = 0; i <= no_of_processes - 1; ++i)
        displacement_vector[i] /= no_of_eq, chunks_per_proc[i] /= no_of_eq;

    MPI_Scatterv(b_matrix, chunks_per_proc, displacement_vector, MPI_DOUBLE, proc_vals, B_SZ, MPI_DOUBLE, root_rank, comm);

    // defining variables
    MPI_Status st;
    REC_SZ = no_of_eq + 1 + 1; // sending no. of eq (n) + pivot (1)+ corresp val of b (1)
    ld recvd_row[REC_SZ];

    ll prev_proc, next_proc, curr, prev_curr;
    curr = rank, prev_curr = -1;
    next_proc = (1 + rank);
    next_proc %= no_of_processes;
    prev_proc = (no_of_processes + rank - 1);
    prev_proc %= no_of_processes;

    // now our task is to get upper triangular matrix and then backsubstitute
    ll piv, cnt = 0;
    int flag = 0;
    MPI_Barrier(comm);
    // fprintf(stderr, "hello %d %lld\n", rank, rows_per_proc);
    // cerr<<endl;
    // cerr<<"hello"<<rank<<rows_per_proc<<endl;
    // we will loop over the rows and communicate in pipelined way
    while (cnt <= rows_per_proc - 1)
    {
        // cout<<endl<<"ko"<<rank<<root_rank<<":"<<cnt<<endl;
        ld send_buf[REC_SZ]; // this will contain first n processed values of rows and next pivot at nth posn and value of b at n+1th index
        // Iterating over all the rows before the current row and
        // previously processed row, to perform elimination
        // corresponding to that row
        ll i = (prev_curr + 1);
        while (i <= curr - 1)
        {
            // cerr<<"ee "<<rank<<cnt<<endl;
            MPI_Recv(recvd_row, REC_SZ, MPI_DOUBLE, prev_proc, i, comm, &st);
            if (curr < (i + no_of_processes - 1))
            {
                flag = 1;
                MPI_Send(recvd_row, REC_SZ, MPI_DOUBLE, next_proc, i, comm);
                // preforming elimination step
            }
            // printf("First");
            Eliminate(i, rank, no_of_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, no_of_processes, permut_variables);
            i++;
        }
        // cerr<<"prev loop"<<rank<<" "<<cnt<<endl;
        // calculating pivot
        ll piv;
        ld mx;
        ll row_id;
        row_id = curr / no_of_processes;
        mx = INT_MIN;
        ll pivot;
        i = curr;
        ld val;
        while (i <= no_of_eq - 1)
        {
            val = proc_rows[i + (no_of_eq * row_id)];
            val = abs(val);
            if (val > mx)
                pivot = i, mx = val;
            i++;
        }
        piv = pivot;
        // ll tem = pivot;
        // printf("First3");
        // performing division step

        swap(proc_rows[pivot + ((curr / no_of_processes) * no_of_eq)], proc_rows[curr + ((curr / no_of_processes) * no_of_eq)]);
        i = curr;
        ld piv_val = proc_rows[curr + ((curr / no_of_processes) * no_of_eq)];
        while (i <= no_of_eq - 1)
        {
            proc_rows[i + ((curr / no_of_processes) * no_of_eq)] /= piv_val;
            i++;
        }
        proc_vals[curr / no_of_processes] /= piv_val;

        // now update the send_buf from proc_rows
        for (ll j = 0; j <= no_of_eq - 1; ++j)
            send_buf[j] = proc_rows[(j + no_of_eq * cnt)];
        send_buf[(no_of_eq + 1)] = proc_vals[cnt], send_buf[no_of_eq] = pivot;
        // but if no_of_processes<2 then last row will also be sent and we dont want that so handling it
        if (no_of_processes >= 2)
        {
            MPI_Send(send_buf, REC_SZ, MPI_DOUBLE, next_proc, curr, comm);
            // cerr<<"dd\n"<<rank<<endl;
        }
        // updating curr and prev_curr
        prev_curr = curr;
        curr = (no_of_processes + curr);
        cnt++;
        Eliminate(prev_curr, rank, no_of_eq, proc_rows, proc_vals, curr, send_buf, rows_per_proc, no_of_processes, permut_variables);
    }
    // recieving non processed rows
    ll i = prev_curr + 1;
    ll lst, prev_lst, count = 1;
    while (i <= no_of_eq - 1)
    {
        MPI_Recv(recvd_row, REC_SZ, MPI_DOUBLE, prev_proc, i, comm, &st);

        if (curr < (i + no_of_processes - 1))
            MPI_Send(recvd_row, REC_SZ, MPI_DOUBLE, next_proc, i, comm);

        Eliminate(i, rank, no_of_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, no_of_processes, permut_variables);

        i++;
    }
    MPI_Barrier(comm);

    // now the back subs part
    lst = no_of_processes * (rows_per_proc)-no_of_processes + rank; // actual last row for curr process
    prev_lst = no_of_eq;
    cnt = rows_per_proc - 1;
    ld res[no_of_eq];
    while (cnt >= 0)
    {
        ld ans;
        for (i = prev_lst - 1; i >= lst + 1; --i)
        {
            ll sd = no_of_eq + i;
            ld x_val;
            MPI_Recv(&x_val, count, MPI_DOUBLE, next_proc, sd, comm, &st);
            if (i < (no_of_processes - 1 + lst))
                MPI_Send(&x_val, count, MPI_DOUBLE, prev_proc, sd, comm);
            for (ll k = cnt; k > -1; --k)
            {
                proc_vals[k] -= (x_val * proc_rows[i + (no_of_eq * k)]);
            }
        }
        ans = proc_vals[cnt];
        for (ll k = cnt - 1; k > -1; --k)
        {
            proc_vals[k] -= (ans * proc_rows[lst + no_of_eq * k]);
        }
        if (no_of_processes >= 2 && lst >= 1)
        {
            MPI_Send(&ans, count, MPI_DOUBLE, prev_proc, no_of_eq + lst, comm);
        }
        prev_lst = lst;
        cnt--;
        lst = lst - no_of_processes;
    }
    MPI_Gatherv(proc_vals, rows_per_proc, MPI_DOUBLE, res, chunks_per_proc, displacement_vector, MPI_DOUBLE, root_rank, comm);

    if (rank == root_rank)
    {
        ll k;
        ld solution[no_of_eq], sol[no_of_eq];
        k = 0;
        for (ll i = 0; i <= no_of_processes - 1; ++i)
        {
            ll j = i;
            while (j <= no_of_eq - 1)
            {
                solution[j] = res[k++];
                j += no_of_processes;
            }
        }
        for (ll i = 0; i <= no_of_eq - 1; ++i)
        {
            ld tv = solution[i];
            sol[permut_variables[i]] = tv;
        }
        FILE *outfile = fopen(argv[2], "w");
        for (ll i = 0; i <= no_of_eq - 1; ++i)
        {
            fprintf(outfile, "%lf ", sol[i]);
        }
        fclose(outfile);
    }

    MPI_Barrier(comm);
    MPI_Finalize();
    return EXIT_SUCCESS;
}