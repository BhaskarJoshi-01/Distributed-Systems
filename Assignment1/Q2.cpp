#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;
int offset_clique_4 = 4;

void combinations(int arr[], int data[], int start, int end, int index, int r, int adjacency_matrix[], int clique_count[], int n)
{
    if (index == r)
    {
        int weight_sum = 0;
        for (int i = 0; i < r; i++)
            for (int j = i + 1; j < r; j++)
                weight_sum += adjacency_matrix[data[i] * n + data[j]];

        if (r == 3 && adjacency_matrix[data[1] * n + data[2]] != -1)
            ++clique_count[weight_sum];
        else if (r == 4)
        {
            if (adjacency_matrix[data[1] * n + data[2]] != -1 &&
                adjacency_matrix[data[0] * n + data[3]] != -1 &&
                adjacency_matrix[data[1] * n + data[3]] != -1 &&
                adjacency_matrix[data[2] * n + data[3]] != -1)
                ++clique_count[offset_clique_4 + weight_sum];
        }
        return;
    }
    for (int i = start; i <= end && end - i + 1 >= r - index; i++)
    {
        data[index] = arr[i];
        combinations(arr, data, i + 1, end, index + 1, r, adjacency_matrix, clique_count, n);
    }
}

void get_cliques_for_vertices(int start_vertex, int end_vertex, int n, int adjacency_matrix[], int clique_count[])
{
    int adj_vert[n];
    int temp[4];
    for (int vertex = start_vertex; vertex <= end_vertex; vertex++)
    {
        int vert_count = 0;
        for (int j = vertex + 1; j < n; j++)
        {
            if (adjacency_matrix[n * vertex + j] != -1)
            {
                adj_vert[vert_count] = j;
                ++vert_count;
            }
        }
        temp[0] = vertex;
        combinations(adj_vert, temp, 0, vert_count - 1, 1, 3, adjacency_matrix, clique_count, n);
        temp[0] = vertex;
        combinations(adj_vert, temp, 0, vert_count - 1, 1, 4, adjacency_matrix, clique_count, n);
    }
}

int main(int argc, char **argv)
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
    if (numprocs >= 12)
    {
        printf("Incorrect Input format \nGive only input as: \n mpirun -np 11 ./a.out <input file> <output file>. \n");
        MPI_Abort(comm, EXIT_FAILURE);
    }

    ifstream fin;
    int n, e;
    if (rank == 0)
    {
        fin.open(argv[1], ios::in);
        fin >> n >> e;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int adjacency_matrix[n * n];

    if (rank == 0)
    {
        memset(adjacency_matrix, -1, sizeof(adjacency_matrix));
        for (int i = 0; i < e; i++)
        {
            int u, v, w;
            fin >> u >> v >> w;
            --u;
            --v;
            adjacency_matrix[u * n + v] = w;
            adjacency_matrix[v * n + u] = w;
        }
        fin.close();
    }
    MPI_Bcast(adjacency_matrix, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    int clique_count[11];
    memset(clique_count, 0, sizeof(clique_count));

    // division on basis of rank
    int size = n / (numprocs - 1);
    int start_vertex = rank * size;
    int end_vertex = min(n - 1, size * (rank + 1) - 1);
    if (rank == numprocs - 1)
        end_vertex = n - 1;
    get_cliques_for_vertices(start_vertex, end_vertex, n, adjacency_matrix, clique_count);

    // gather and reduce
    int result[11];
    MPI_Reduce(clique_count, result, 11, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        int k = 0;
        freopen(argv[2], "w", stdout);

        for (int i = 0; i < 4; i++)
            cout << "3 " << i << " " << result[k++] << endl;

        for (int i = 0; i < 7; i++)
            cout << "4 " << i << " " << result[k++] << endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS; // defined as 0 in stdlib
}