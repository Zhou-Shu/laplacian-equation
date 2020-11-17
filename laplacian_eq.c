#include "laplacian_eq.h"

int main(int argc, char *argv[])
{
    double *T;
    int size, rank, opt;
    int period = 0;
    double delta;
    struct timeval start, end;
    MPI_Comm link_comm;

    //  initialize Temperature matrix
    T = init_T(top_value, bottom_value, left_value, right_value);

    while ((opt = getopt_long(argc, argv, optstr, opts, NULL)) != -1)
    {
        switch (opt)
        {
        case 'p': // parallel part
            gettimeofday(&start, NULL);

            MPI_Init(&argc, &argv);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            // set a reasonable physical topology for this task
            MPI_Cart_create(MPI_COMM_WORLD, 1, &size, &period, 1, &link_comm);
            solve_laplacian_eq_parallel(T, link_comm);
            MPI_Finalize();

            gettimeofday(&end, NULL);
            delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            if (rank == 0)
            {
                if (optarg != NULL)
                    print_array(T);
                printf("time=%12.10f\n", delta);
            }
            return 0;
            break;
        case 's': // serial part
            gettimeofday(&start, NULL);

            solve_laplacian_eq_serial(T);

            gettimeofday(&end, NULL);
            delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            if (optarg != NULL)
                print_array(T);
            printf("time=%12.10f\n", delta);
            return 0;
            break;
        default:
            printf("option requires: -p/-s/--parallel/--serial\n");
            return 1;
            break;
        }
    }
    printf("option requires: -p/-s/--parallel/--serial\n");
        return 1;
}

// initialize the matrix, set the boundary condition

double *init_T(int top, int bottom, int left, int right)
{
    double *T = (double *)calloc(rows_num * cols_num, sizeof(double));

    for (size_t i = 0; i < cols_num; i++)
    {
        T[i] = top;
        T[(rows_num - 1) * cols_num + i] = bottom;
    }
    for (size_t i = 0; i < rows_num; i++)
    {
        T[i * cols_num] = left;
        T[(i)*cols_num + cols_num - 1] = right;
    }
    return T;
}

// using serial algorithm to solve the eq

void solve_laplacian_eq_serial(double *T)
{
    for (size_t i = 0; i < iterate; i++)
    {
        for (size_t x = 1; x < rows_num - 1; x++)
        {
            for (size_t y = 1; y < cols_num - 1; y++)
            {
                T[x * cols_num + y] = (T[(x + 1) * cols_num + y] + T[(x - 1) * cols_num + y] + T[x * cols_num + y + 1] + T[x * cols_num + y - 1]) * 0.25;
            }
        }
    }
}

// using parallel algorithm to solve the eq

void solve_laplacian_eq_parallel(double *T, MPI_Comm link_comm)
{
    int right_nbr, left_nbr, rank, size;
    int top_tag = 1;
    int bottom_tag = 2;
    MPI_Status status;

    MPI_Comm_rank(link_comm, &rank);
    MPI_Comm_size(link_comm, &size);
    MPI_Cart_shift(link_comm, 0, 1, &left_nbr, &right_nbr);

    int *index_array = (int *)calloc(size + 1, sizeof(int));
    int *size_array = (int *)calloc(size + 1, sizeof(int));
    int *disp_array = (int *)calloc(size + 1, sizeof(int));
    split_matrix(size, index_array, size_array, disp_array);

    int start_index = index_array[rank];
    int end_index = index_array[rank + 1];

    for (size_t i = 0; i < iterate; i++)
    {

        if (rank == 0)
        {
            for (size_t x = start_index + 1; x < end_index; x++)
            {
                for (size_t y = 1; y < cols_num - 1; y++)
                {
                    T[x * cols_num + y] = (T[(x + 1) * cols_num + y] + T[(x - 1) * cols_num + y] + T[x * cols_num + y + 1] + T[x * cols_num + y - 1]) * 0.25;
                }
            }
            MPI_Send(&T[(end_index - 1) * cols_num], cols_num, MPI_DOUBLE, right_nbr, bottom_tag, link_comm);
            MPI_Recv(&T[end_index * cols_num], cols_num, MPI_DOUBLE, right_nbr, top_tag, link_comm, &status);
        }

        // processes with even rank send msg first to avoid the deadlock

        else if (rank % 2 == 0)
        {
            for (size_t x = start_index; x < end_index; x++)
            {
                for (size_t y = 1; y < cols_num - 1; y++)
                {
                    T[x * cols_num + y] = (T[(x + 1) * cols_num + y] + T[(x - 1) * cols_num + y] + T[x * cols_num + y + 1] + T[x * cols_num + y - 1]) * 0.25;
                }
            }
            MPI_Send(&T[start_index * cols_num], cols_num, MPI_DOUBLE, left_nbr, top_tag, link_comm);
            if (rank != size - 1)
                MPI_Send(&T[(end_index - 1) * cols_num], cols_num, MPI_DOUBLE, right_nbr, bottom_tag, link_comm);
            MPI_Recv(&T[(start_index - 1) * cols_num], cols_num, MPI_DOUBLE, left_nbr, bottom_tag, link_comm, &status);
            if (rank != size - 1)
                MPI_Recv(&T[end_index * cols_num], cols_num, MPI_DOUBLE, right_nbr, top_tag, link_comm, &status);
        }

        //Also, processes with odd rank recv msg first

        else if (rank % 2 == 1)
        {
            if (rank != size - 1)
                MPI_Recv(&T[end_index * cols_num], cols_num, MPI_DOUBLE, right_nbr, top_tag, link_comm, &status);
            MPI_Recv(&T[(start_index - 1) * cols_num], cols_num, MPI_DOUBLE, left_nbr, bottom_tag, link_comm, &status);
            if (rank != size - 1)
                MPI_Send(&T[(end_index - 1) * cols_num], cols_num, MPI_DOUBLE, right_nbr, bottom_tag, link_comm);
            MPI_Send(&T[start_index * cols_num], cols_num, MPI_DOUBLE, left_nbr, top_tag, link_comm);

            for (size_t x = start_index; x < end_index; x++)
            {
                for (size_t y = 1; y < cols_num - 1; y++)
                {
                    T[x * cols_num + y] = (T[(x + 1) * cols_num + y] + T[(x - 1) * cols_num + y] + T[x * cols_num + y + 1] + T[x * cols_num + y - 1]) * 0.25;
                }
            }
        }
    }

    //gather all the informations to the master

    // printf("rank=%d, start from %d, end to %d, will send data from %d, send %d data\n", rank,index_array[rank] , index_array[rank+1],disp_array[rank], size_array[rank]);
    MPI_Gatherv(&T[start_index * cols_num], cols_num * (end_index - start_index), MPI_DOUBLE, T, size_array, disp_array, MPI_DOUBLE, 0, link_comm);
}

// try to use a more reasonable way to split the task for every processes

void split_matrix(int num_process, int *index_array, int *size_array, int *disp_array)
{
    for (size_t i = 0; i < rows_num; i++)
    {

        index_array[i % num_process + 1]++;
        size_array[i % num_process] += cols_num;
        disp_array[i % num_process + 1] += cols_num;
    }
    for (size_t i = 1; i < num_process + 1; i++)
    {
        index_array[i] += index_array[i - 1];
        disp_array[i] += disp_array[i - 1];
    }
    index_array[num_process]--;
}

// print the result

void print_array(double *T)
{

    for (size_t i = 0; i < rows_num; i++)
    {
        for (size_t j = 0; j < cols_num; j++)
        {
            printf("%.1f ", T[i * cols_num + j]);
        }
        printf("\n");
    }
}
