/* Include benchmark-specific header. */
#include "jacobi-2d.h"
#include <mpi.h>
double bench_t_start, bench_t_end;

static
void init_array (int n,
   double A[n][n],
   double B[n][n])
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
 A[i][j] = ((double) i*(j+2) + 2) / n;
 B[i][j] = ((double) i*(j+3) + 3) / n;
      }
}

static
void kernel_jacobi_2d(int tsteps,
       int n,
       double A[ n][n],
       double B[ n][n])
{
  int t, i, j, k, start, fin;
  double ibuffer[n-2], obuffer[n-2];
  int np, myid;
  MPI_Comm world;
  MPI_Status s;
  
  
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world, &np);
  MPI_Comm_rank(world, &myid);
  
  k = (n-2) / np;
  if ((n-2) % np) k++; 
  if (k == 0) k = 1;
  start = 1+myid*k;
  fin = start + k;
  if (start >= n-1) fin = 0;
  if (fin > n-1) fin = n-1;
  
  printf("I'm %d and my start = %d; fin = %d\n", myid, start, fin);
  
  for (t = 0; t < tsteps; t++)
    {
      for (i = start; i < fin; i++)
        for (j = 1; j < n - 1; j++)
            B[i][j] = 0.2 * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
      
      if (fin)
      {
        for (i = 1; i < n-1; i++) obuffer[i] = B[start][i];
        if (start > 1)
            MPI_Send ( &obuffer, n-2, MPI_DOUBLE, myid-1, 1, world);
        if (fin < n-1)
        {
            MPI_Recv ( &ibuffer, n-2, MPI_DOUBLE, myid+1, 1, world, &s);
            for (i = 1; i < n-1; i++) B[fin][i] = ibuffer[i];
        }
        for (i = 1; i < n-1; i++) obuffer[i] = B[fin-1][i];
        if (fin < n-1)
            MPI_Send ( &obuffer, n-2, MPI_DOUBLE, myid+1, 2, world);
        if (start > 1)
        {
            MPI_Recv ( &ibuffer, n-2, MPI_DOUBLE, myid-1, 2, world, &s);
            for (i = 1; i < n-1; i++) B[start-1][i] = ibuffer[i];
        }
      }
      MPI_Barrier(world);
      
      for (i = start; i < fin; i++)
        for (j = 1; j < n - 1; j++)
            A[i][j] = 0.2 * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
      
      if (fin)
      {
        for (i = 1; i < n-1; i++) obuffer[i] = A[start][i];
        if (start > 1)
            MPI_Send ( &obuffer, n-2, MPI_DOUBLE, myid-1, 1, world);
        if (fin < n-1)
        {
            MPI_Recv ( &ibuffer, n-2, MPI_DOUBLE, myid+1, 1, world, &s);
            for (i = 1; i < n-1; i++) A[fin][i] = ibuffer[i];
        }
        for (i = 1; i < n-1; i++) obuffer[i] = A[fin-1][i];
        if (fin < n-1)
            MPI_Send ( &obuffer, n-2, MPI_DOUBLE, myid+1, 2, world);
        if (start > 1)
        {
            MPI_Recv ( &ibuffer, n-2, MPI_DOUBLE, myid-1, 2, world, &s);
            for (i = 1; i < n-1; i++) A[start-1][i] = ibuffer[i];
        }
      }
      MPI_Barrier(world);
      
    }
}


int main(int argc, char** argv)
{
  int n = N;
  int tsteps = TSTEPS;
  double (*A)[n][n]; A = (double(*)[n][n])malloc ((n) * (n) * sizeof(double));
  double (*B)[n][n]; B = (double(*)[n][n])malloc ((n) * (n) * sizeof(double));

  init_array (n, *A, *B);

  MPI_Init(&argc, &argv);

  bench_t_start = MPI_Wtime();

  kernel_jacobi_2d(tsteps, n, *A, *B);

  bench_t_end = MPI_Wtime();
  printf ("Time in seconds = %f\n", bench_t_end - bench_t_start);
  
  free((void*)A);
  free((void*)B);
  
  MPI_Finalize();
  return 0;
}
