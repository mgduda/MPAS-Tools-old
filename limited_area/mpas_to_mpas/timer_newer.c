#include <stdlib.h>
#include <stdio.h>
#include <time.h>

struct timespec start_time[10];
struct timespec end_time[10];

void start_timer(int n)
{
   clock_gettime(CLOCK_MONOTONIC_RAW, &start_time[n]);
}

void stop_timer(int n, int * secs, int * n_secs)
{
   clock_gettime(CLOCK_MONOTONIC_RAW, &end_time[n]);

   *secs = (int)(end_time[n].tv_sec - start_time[n].tv_sec);
   *n_secs = (int)(end_time[n].tv_nsec - start_time[n].tv_nsec);

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
}

