//#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

int myrank;

void *
myalloc(size_t size)
{
  size_t align = 64;
#if 0
  char *a = malloc(size+32);
  int o = ((size_t)a) & 31;
  o = (32 - o) & 31;
  char *b = a + o;
  return (void *)b;
#else
  void *b = NULL;
  int err = posix_memalign(&b, align, size);
  if(err) {
    printf("posix_memalign(&%p, %li, %li) failed!\n", b, align, size);
    exit(1);
  }
  return b;
#endif
}

void
threadSum(int n, ThreadCtx *TC)
{
  int tid = TC->tid;
  int nid = TC->nid;
  int volatile *tbar = TC->tbar;
  int tb1 = tbar[tid] + 1;
  int b = 1;
  while((tid&b)==0) {
    int t1 = tid + b;
    if(t1>=nid) break;
    while(tbar[t1]<tb1);
    SYNC;
    double *x = TC->sum[tid];
    double *y = TC->sum[t1];
    for(int i=0; i<n; i++) x[i] += y[i];
    b *= 2;
  }
  SYNC;
  tbar[tid] = tb1;
}

void
gsum(double *x, int n, ThreadCtx *TC)
{
  int tid = TC->tid;
  //printf("tid: %i  x: %g\n", tid, *x);
  TC->sum[tid] = x;
  threadSum(n, TC);
  if(tid==0) {
    QMP_sum_double_array(x, n);
  }
  TWAIT0;
  if(tid!=0) {
    double *y = TC->sum[0];
    for(int i=0; i<n; i++) x[i] = y[i];
  }
  T0WAIT;
  //printf("tid: %i  x: %g\n", tid, *x);
  //printf0("gsum:");
  //for(int i=0; i<n; i++) printf0(" %g", x[i]);
  //printf0("\n");
}
