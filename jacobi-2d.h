#ifndef _JACOBI_2D_H
#define _JACOBI_2D_H 
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#define MINI_DATASET
# endif
# if !defined(TSTEPS) && !defined(N)
# ifdef MINI_DATASET
#define TSTEPS 20
#define N 30
# endif
# ifdef SMALL_DATASET
#define TSTEPS 40
#define N 90
# endif
# ifdef MEDIUM_DATASET
#define TSTEPS 100
#define N 250
# endif
# ifdef LARGE_DATASET
#define TSTEPS 500
#define N 1300
# endif
# ifdef EXTRALARGE_DATASET
#define TSTEPS 1000
#define N 2800
# endif
#endif
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#endif
