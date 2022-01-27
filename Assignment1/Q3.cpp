#include <mpi.h>
#include <bits/stdc++.h>
#include <fstream>
#define ll long long int
#define ld double
#define fw(i, s, e) for (ll i = s; i < e; ++i) // forward loop so fw
#define fe(i, s, e) for (ll i = s; i <= e; ++i)
#define fb(i, e, s) for (ll i = e; i >= s; --i) // for loop but backward so fb
using namespace std;
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

void perform_elimination(int recv_id, int id, int num_eq, ld *proc_rows, ld *proc_vals, int curr, ld *recvd_row, int rows_per_proc, int num_proc, int *var_perm)
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
                ld piv_val = proc_rows[(num_eq * i) + recv_id];
                for (ll j = 0; j <= num_eq - 1; ++j)
                    proc_rows[j + (num_eq * i)] -= (recvd_row[j] * piv_val);
                proc_vals[i] = proc_vals[i] - recvd_row[num_eq + 1] * piv_val;
            }
        }

        i++;
    }
}


ll compute_pivot(int curr, int num_proc, int num_eq, ld *proc_rows)
{
    ld mx;
    ll row_id;
    row_id = curr / num_proc;
    mx = INT_MIN;
    ll pivot, i = curr;
    ld val;
    while (i <= num_eq - 1)
    {
        val = proc_rows[i + (num_eq * row_id)];
        val = abs(val);
        if (val > mx)
            pivot = i, mx = val;
        i++;
    }
    return pivot;
}

void perform_division(int id, int curr, ld *proc_rows, int pivot, int num_proc, int num_eq, int rows_per_proc, ld *proc_vals)
{

    swap(proc_rows[pivot + ((curr / num_proc) * num_eq)], proc_rows[curr + ((curr / num_proc) * num_eq)]);
    ll i = curr;
    ld piv_val = proc_rows[curr + ((curr / num_proc) * num_eq)];
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
    ll num_eq=0;
    // reading number of equations from file first and creating
    // variables accordingly for taking matrix as input and distributing

    if (id == root_rank)
    {
        inputfp = fopen(argv[1], "r");
        fscanf(inputfp, "%lld", &num_eq);
    }
    printf("First4");

    ll oned_size = pow(num_eq, 2);
    // converting n*n A matrix that was input to 1*n matrix

    ld eq_mat[oned_size], val_mat[num_eq];
    // sending no. of equations to every process
    MPI_Bcast(&num_eq, 1, MPI_INT, root_rank, comm);
    printf("First5 %d\n", id);

    // Defining variables for dividing processes
    ll rpp, rows_per_proc, np0, np1;
    int var_perm[num_eq]; // stores the column no. in array
    for (ll i = 0; i < num_eq; ++i)
        var_perm[i] = i;
    rpp = num_eq / num_proc; // no of process that shoud be divided into processes
    // like 12/4=3

    // now basically what i m trying to do is pair conditions like
    // 15/4 so there will be like 3 processes with one more row and only
    // one process with 15/4

    ll A_SZ, B_SZ, REC_SZ, rem_proc = num_eq % num_proc;
    np0 = num_proc - rem_proc;
    np1 = num_proc - np0;
    // if process rank is less than np1 so we should give it one more pocess
    rows_per_proc = ((id < np1) ? (1 + rpp) : rpp);
    fprintf(stderr, "hello %d %lld\n", id, rows_per_proc);
    int divs[num_proc], displs[num_proc];
    memset(divs, 0, sizeof(divs));
    if (id == root_rank)
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
            cerr<<"\ndivs update"<<divs[i % num_proc]<<endl;
            // now we have assigned row in process , now we just read and allocate from imput
            while (k < num_eq)
            {
                fscanf(inputfp, "%lf", &eq_mat[k + num_eq * eff_row]);
                k++;
            }
            fscanf(inputfp, "%lf", &val_mat[eff_row]);
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
    // next we have distributed the A matrix to other porcesses
    A_SZ = rows_per_proc * num_eq;
    B_SZ = A_SZ / num_eq;
    ld proc_rows[A_SZ]; // part of A matrix for a process
    ld proc_vals[B_SZ]; // part of B matrix for a process

    // MPI_Scatterv is in which the data dispatched from the root process
    // can vary in the number of elements, and the location from which load
    // these elements in the root process buffer.
    if(id==0){

    // for (int i = 0; i < A_SZ; i++)
    // {
    //     cerr << eq_mat[i] << " ";
    // }
    // cerr << "\n";
    cerr<<"\n"<<id<<"DIVS: ";
    for (int i = 0; i < num_proc; i++)
    {
        cerr << divs[i] << " ";
    }
    cerr<<endl;
    cerr << "disps\n";
    for (int i = 0; i < num_proc; i++)
    {
        cerr << displs[i] << " ";
    }
    cerr << "\n";
    }
    MPI_Scatterv(eq_mat, divs, displs, MPI_DOUBLE, proc_rows, A_SZ, MPI_DOUBLE, root_rank, comm);
    cerr << "hi " << id << endl;

    // similarly we need to do for B matrix
    // as divs and disps were multiplied by num_eq above
    // changing them back
    for (ll i = 0; i <= num_proc - 1; ++i)
        displs[i] /= num_eq, divs[i] /= num_eq;

    MPI_Scatterv(val_mat, divs, displs, MPI_DOUBLE, proc_vals, B_SZ, MPI_DOUBLE, root_rank, comm);

    // defining variables
    MPI_Status st;
    REC_SZ = num_eq + 1 + 1; // sending no. of eq (n) + pivot (1)+ corresp val of b (1)
    ld recvd_row[REC_SZ];

    ll prev_proc, next_proc, curr, prev_curr;
    curr = id, prev_curr = -1;
    next_proc = (1 + id);
    next_proc %= num_proc;
    prev_proc = (num_proc + id - 1);
    prev_proc %= num_proc;

    // now our task is to get upper triangular matrix and then backsubstitute
    ll piv, cnt = 0;
    int flag = 0;
    MPI_Barrier(comm);
    fprintf(stderr, "hello %d %lld\n", id, rows_per_proc);
    cerr<<endl;
    // cerr<<"hello"<<id<<rows_per_proc<<endl;
    // we will loop over the rows and communicate in pipelined way
    while (cnt <= rows_per_proc - 1)
    {
    cout<<endl<<"ko"<<id<<root_rank<<":"<<cnt<<endl;
        ld send_buf[REC_SZ]; // this will contain first n processed values of rows and next pivot at nth posn and value of b at n+1th index
        // Iterating over all the rows before the current row and
        // previously processed row, to perform elimination
        // corresponding to that row
        ll i = prev_curr + 1;
        while (i <= curr - 1)
        {
            cerr<<"ee "<<id<<cnt<<endl;
            MPI_Recv(recvd_row, REC_SZ, MPI_DOUBLE, prev_proc, i, comm, &st);
            if (curr < (i + num_proc - 1))
                MPI_Send(recvd_row, REC_SZ, MPI_DOUBLE, next_proc, i, comm);
            // preforming elimination step
            flag = 1;
            printf("First");
            perform_elimination(i, id, num_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, num_proc, var_perm);
            i++;
        }
        cerr<<"prev loop"<<id<<" "<<cnt<<endl;

        piv = compute_pivot(curr, num_proc, num_eq, proc_rows);
        ll tem = piv;
        printf("First3");

        perform_division(id, curr, proc_rows, tem, num_proc, num_eq, rows_per_proc, proc_vals);
        // now update the send_buf from proc_rows
        for (ll j = 0; j <= num_eq - 1; ++j)
            send_buf[j] = proc_rows[(j + num_eq * cnt)];
        send_buf[(num_eq + 1)] = proc_vals[cnt], send_buf[num_eq] = tem;
        // but if num_proc<2 then last row will also be sent and we dont want that so handling it
        if (num_proc >= 2){
            MPI_Send(send_buf, REC_SZ, MPI_DOUBLE, next_proc, curr, comm);
            cerr<<"dd\n"<<id<<endl;
        }
        // updating curr and prev_curr
        prev_curr = curr;
        curr = (num_proc + curr);
        cnt++;
        perform_elimination(prev_curr, id, num_eq, proc_rows, proc_vals, curr, send_buf, rows_per_proc, num_proc, var_perm);
    }
    // recieving non processed rows
    ll i = prev_curr + 1;
    ll lst, prev_lst, count = 1;
    while (i <= num_eq - 1)
    {
        MPI_Recv(recvd_row, REC_SZ, MPI_DOUBLE, prev_proc, i, comm, &st);

        if (curr < (i + num_proc - 1))
            MPI_Send(recvd_row, REC_SZ, MPI_DOUBLE, next_proc, i, comm);

        perform_elimination(i, id, num_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, num_proc, var_perm);

        i++;
    }
    MPI_Barrier(comm);

    // now the back subs part
    lst = num_proc * (rows_per_proc)-num_proc + id; // actual last row for curr process
    prev_lst = num_eq;
    cnt = rows_per_proc - 1;
    ld res[num_eq];
    while (cnt >= 0)
    {
        ld ans;
        for (i = prev_lst - 1; i >= lst + 1; --i)
        {
            ll sd = num_eq + i;
            ld x_val;
            MPI_Recv(&x_val, count, MPI_DOUBLE, next_proc, sd, comm, &st);
            if (i < (num_proc - 1 + lst))
                MPI_Send(&x_val, count, MPI_DOUBLE, prev_proc, sd, comm);
            for (ll k = cnt; k > -1; --k)
            {
                proc_vals[k] -= (x_val * proc_rows[i + (num_eq * k)]);
            }
        }
        ans = proc_vals[cnt];
        for (ll k = cnt - 1; k > -1; --k)
        {
            proc_vals[k] -= (ans * proc_rows[lst + num_eq * k]);
        }
        if (num_proc >= 2 && lst >= 1)
        {
            MPI_Send(&ans, count, MPI_DOUBLE, prev_proc, num_eq + lst, comm);
        }
        prev_lst = lst;
        cnt--;
        lst = lst - num_proc;
    }
    MPI_Gatherv(proc_vals, rows_per_proc, MPI_DOUBLE, res, divs, displs, MPI_DOUBLE, root_rank, comm);

    if (id == root_rank)
    {
        ll k;
        ld solution[num_eq], sol[num_eq];
        k = 0;
        for (ll i = 0; i <= num_proc - 1; ++i)
        {
            ll j = i;
            while (j <= num_eq - 1)
            {
                solution[j] = res[k++];
                j += num_proc;
            }
        }
        for (ll i = 0; i <= num_eq - 1; ++i)
        {
            ld tv = solution[i];
            sol[var_perm[i]] = tv;
        }
        FILE *outfile = fopen(argv[2], "w");
        for (ll i = 0; i <= num_eq - 1; ++i)
        {
            fprintf(outfile, "%lf ", sol[i]);
        }
        fclose(outfile);
    }

    MPI_Barrier(comm);
    MPI_Finalize();
    return EXIT_SUCCESS;
}