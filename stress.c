#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "common.h"

static Layout layout;

double
dtime(void)
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
#if 0
  struct timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return tp.tv_sec + 1e-9*tp.tv_nsec;
#else
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec + 1e-6*tp.tv_usec;
#endif
#endif
}

void
init(int *argc, char ***argv)
{
  QMP_thread_level_t provided;
  int err = QMP_init_msg_passing(argc, argv, QMP_THREAD_FUNNELED, &provided);
  if(err!=0) {
    printf("%s\n", QMP_error_string(err));
    printf("com_qmp: Initialize QMP failed.\n");
    fflush(stdout);
    exit(err);
  }
}

void
fini(void)
{
  QMP_finalize_msg_passing();
}

void
dump(real *x, int n, int vl)
{
  for(int i=0; i<n; i++) {
    printf("%12g", x[i*vl]);
    for(int j=1; j<vl; j++) {
      printf(" %12g", x[i*vl+j]);
    }
    printf("\n");
  }
}

void
scalels(int *ls, int nd, int n)
{
  while(n>1) {
    for(int f=2; f<=n; f++) {
      if((n/f)*f==n) {
	int k = 0;
	for(int i=1; i<nd; i++) if(ls[i]<=ls[k]) k = i;
	ls[k] *= f;
	n /= f;
      }
    }
  }
}

int
main(int argc, char *argv[])
{
  double t0 = dtime();
  init(&argc, &argv);
  t0 = dtime() - t0;

#define ND 4
  int nd = ND;
  int nd2 = 2*nd;
  int ls[ND];
  double nrepfac = 2e7;
  for(int i=0; i<nd; i++) ls[i] = 12;
  for(int i=1,j=0; argv[i]; i++) {
    if(argv[i][0]=='r') nrepfac = atof(argv[i]+1);
    else
      ls[j++] = atoi(argv[i]);
  }

  layout.nranks = QMP_get_number_of_nodes();
  layout.myrank = QMP_get_node_number();
  scalels(ls, nd, layout.nranks);
  myrank = layout.myrank;
  printf0("init time: %g\n", t0);
  t0 = dtime();
  P(stagDslashSetup)(&layout, nd, ls, NULL);
  t0 = dtime() - t0;
  printf0("stagDslashSetup time: %g\n", t0);
  int nthreads = 1;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }

  int vol = layout.nSites;
  int nc = 3;
  int nelemV = 2*nc;
  int nelemM = 2*nc*nc;
  int veclen = layout.nSitesInner;
  Field v1, v2, v3, u[nd2];
  Subset all;
  layoutSubset(&all, &layout, "all");

  fieldNew(&v1, &layout, nelemV);
  fieldNew(&v2, &layout, nelemV);
  fieldNew(&v3, &layout, nelemV);
  for(int i=0; i<nd2; i++) {
    fieldNew(&u[i], &layout, nelemM);
  }

  double nv1;
  beginParallelSubset(all) {
    F_eq_r(&v1, 0);
    for(int i=ti0; i<ti1; i++) {
      int k = i*nelemV*veclen;
      for(int j=0; j<nelemV*veclen; j++) {
	v1.f[k+j] = 0*k + (j/veclen)+1;
      }
    }
    F_eq_r(&v2, 0);
    F_eq_r(&v3, 0);
    for(int d=0; d<nd2; d++) {
      F_eq_r(&u[d], 0);
      for(int i=ti0; i<ti1; i++) {
	int k = i*nelemM*veclen;
	for(int j=0; j<nelemM*veclen; j++) {
	  u[d].f[k+j] = 0*k + (j/veclen)+1;
	}
#if 0
	for(int j1=0; j1<nc; j1++) {
	  for(int j2=0; j2<nc; j2++) {
	    if(j1>=j2) {
	      int j0 = k + (j1+j2*nc)*2*veclen;
	      for(int l=0; l<veclen; l++) {
		u[d].f[j0+l] = 1+j1;
	      }
	    }
	  }
	}
#endif
      }
      P(makeU)(u[d].f, nc, ti0, ti1);
    }
    double r;
    r = lnorm2_F(&v1);
    sumDoubleArray(&r, 1);
    singleThread { nv1 = r; }
  } endParallelSubset;
  printf0("%-20.16g\n", nv1);

  X_eq_r(v2.f, 0, 2*nc, 0, 1);
  P(X_peq_M_times_xXp)(v2.f, u[0].f, (int[1]){0}, 0, v1.f, 2*nc, nc, 2, 0, 1);
  double check = norm2_X(v2.f, 2*nc, 0, 1);
  //double check = norm2_X(u[0].f, 2*nc*nc, 0, 1);
  printf0("%-20.16g\n", check);
  //dump(v1.f, nelemV, veclen);
  //dump(u[0].f, nelemM, veclen);
  //dump(v2.f, nelemV, veclen);
  //X_eq_r(v2.f, 0, 2*nc, 0, 2);
  //P(X_peq_M_times_xXp)(v2.f, u[0].f, (int[3]){0,1,2}, 0, v1.f, 2*nc, nc, 2, 0, 2);
  //check = norm2_X(v2.f, 2*nc, 0, 2);
  //check = norm2_X(u[1].f, 2*nc*nc, 0, 2);
  //printf("%-20.16g\n", check);
  //dump(v1.f, 2*nelemV, veclen);
  //dump(u[0].f, 2*nelemM, veclen);
  //dump(v2.f, 2*nelemV, veclen);

  Transporter T[nd2];
  for(int d=0; d<nd2; d++) {
    ShiftIndices *si = P(getShift)(d/2, 1-2*(d&1), "all");
    //transporterNew(&T0, NULL, si0, &all, nelemV, nc);
    transporterNew(&T[d], u[d].f, si, &all, nelemV, nc);
  }

  int nrep = (int)(1 + (nrepfac*nthreads)/vol);
  nrep = 10000;
  printf0("nrep: %i\n", nrep);

  for(int its=0; its<3600; its++) {
    beginParallelSubset(all) {
      F_eq_r(&v2, 0);
    } endParallelSubset;
    t0 = dtime();
    for(int i=0; i<nrep; i++) {
      int ns = 8;
      beginParallelSubset(all) {
	for(int d=0; d<ns; d++) {
	  transporterStartT(&T[d], v2.f, v1.f, 2, TC);
	}
	for(int d=0; d<ns; d++) {
	  transporterLocalT(&T[d], TC);
	}
	for(int d=0; d<ns; d++) {
	  transporterWaitT(&T[d], TC);
	}
      } endParallelSubset;
    }
    t0 = dtime() - t0;
    double flops = 8*8*nc*nc*vol;
    printf0("xportV: %10g %10g mf: %g\n",t0,1e6*t0/nrep,1e-6*flops*nrep/t0);
    beginParallelSubset(all) {
      double r;
      r = lnorm2_F(&v2);
      sumDoubleArray(&r, 1);
      singleThread { nv1 = r; }
    } endParallelSubset;
    double r = 0, tol = 0;
    if(its==0) {
      r = nv1/check/layout.nSitesOuter/nrep/nrep/(nd2*nd2)/layout.nranks;
      check = nv1;
      tol = 1e-10;
    } else {
      r = nv1/check;
      tol = 1e-20;
    }
    printf0("%-20.16g\t%-20.16g\n", nv1, r);
    if(fabs(r-1)>tol) {
      printf("failure on rank %i\n", layout.myrank);
      QMP_abort(-1);
    }
  }

  fini();
  return 0;
}
