#ifndef _LAPLACIAN_EQ_
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

#define rows_num 1000   // # of rows of the Temperature matrix
#define cols_num 1000  // # of cols of the Temperature matrix
#define iterate 1000    // # of iteration
#define top_value 20    // the boundary condition of top line
#define bottom_value 40 // the boundary condition of bottom line
#define left_value 40   // the boundary condition of left line
#define right_value 40  // the boundary condition of right line

double *init_T(int, int, int, int);
void solve_laplacian_eq_serial(double *);
void solve_laplacian_eq_parallel(double *, MPI_Comm);
void print_array(double *);
void split_matrix(int,int* ,int*,int*);

// command line's options

char *optstr = "p::s::";
struct option opts[] = {
    {"parallel", optional_argument, NULL, 'p'},
    {"serial", optional_argument, NULL, 's'},
    {0, 0, 0, 0},
};

// serial execution: mpiexec -n 1 laplacian_eq.exe --serial / -s
// parallel execution: mpiexec -n 5 laplacian_eq.exe --parallel /-p
// serial execution and print results: mpiexec -n 1 laplacian_eq.exe --serial=p /-sp
// parallel execution and print results:: mpiexec -n 5 laplacian_eq.exe --parallel=p /-pp

#endif
