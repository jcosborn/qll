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

int
main(int argc, char *argv[])
{
  double t0 = dtime();
  qopStagDslashInit(&argc, &argv);
  t0 = dtime() - t0;

#define ND 4
  int nd = ND;
  int ls[ND];
  double nrepfac = 1e8;
  for(int i=0; i<nd; i++) ls[i] = 4;
  for(int i=1,j=0; argv[i]; i++) {
    if(argv[i][0]=='r') nrepfac = atof(argv[i]+1);
    else
      ls[j++] = atoi(argv[i]);
  }

  layout.nranks = QMP_get_number_of_nodes();
  layout.myrank = QMP_get_node_number();
  myrank = layout.myrank;
  printf0("qopStagDslashInit time: %g\n", t0);
  t0 = dtime();
  P(stagDslashSetup)(&layout, nd, ls, NULL);
  t0 = dtime() - t0;
  printf0("stagDslashSetup time: %g\n", t0);
  qopStagDslashSetup(&layout, NULL);
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
  Field v1, v2, v3, u, ua, u1, u2, u3, u4;
  Subset all;
  layoutSubset(&all, &layout, "all");

  fieldNew(&v1, &layout, nelemV);
  fieldNew(&v2, &layout, nelemV);
  fieldNew(&v3, &layout, nelemV);
  fieldNew(&u, &layout, nelemM);
  fieldNew(&ua, &layout, nelemM);
  fieldNew(&u1, &layout, nelemM);
  fieldNew(&u2, &layout, nelemM);
  fieldNew(&u3, &layout, nelemM);
  fieldNew(&u4, &layout, nelemM);

  double nv1, nu1;
  beginParallelSubset(all) {
    F_eq_r(&v1, 0);
    for(int i=ti0; i<ti1; i++) {
      int k = i*nelemV*veclen;
      for(int j=0; j<nelemV*veclen; j++) {
	v1.f[k+j] = k + j;
      }
    }
    F_eq_r(&v2, 0);
    F_eq_r(&u, 0);
    F_eq_r(&u1, 0);
    for(int i=ti0; i<ti1; i++) {
      int k = i*nelemM*veclen;
      for(int j=0; j<nelemM*veclen; j++) {
	u1.f[k+j] = k + j;
      }
      for(int j1=0; j1<nc; j1++) {
	for(int j2=0; j2<nc; j2++) {
	  //if(j1==j2) {
	  int j0 = k + (j1+j2*nc)*2*veclen;
	  for(int l=0; l<veclen; l++) {
	    //u.f[j0+l] = 1;
	    //u.f[j0+l] = 1-(j1*j1)*(j2&1)+(j1&2)*(j2&2);
	    //u.f[j0+veclen+l] = 1+j1+j2;
	    u.f[j0+l] = 1+j0/(double)(1+j1+j2);
	    u.f[j0+veclen+l] = 1+j0/(double)(1+j1*j2);
	  }
	  //}
	}
      }
    }
    P(makeU)(u.f, nc, ti0, ti1);
    F_eq_r(&u2, 0);
    double r[2];
    r[0] = lnorm2_F(&v1);
    r[1] = lnorm2_F(&u1);
    sumDoubleArray(r, 2);
    singleThread { nv1 = r[0]; nu1 = r[1]; }
  } endParallelSubset;
  int dir = 3;
  gaugeDoublestore((Field*[1]){&ua}, (Field*[1]){&u},
		   (int[1]){dir}, (int[1]){-1}, 1);
  //dump(u.f, nelemM, veclen);
  //dump(ua.f, nelemM, veclen);

  Transporter T0, M0, Tb0;
  ShiftIndices *si0 = P(getShift)(dir, 1, "all"); // from forward 0
  ShiftIndices *sib0 = P(getShift)(dir, -1, "all"); // from forward 0
  layoutSubset(&all, &layout, "all");
  //transporterNew(&T0, NULL, si0, &all, nelemV, nc);
  //transporterNew(&Tb0, NULL, sib0, &all, nelemV, nc);
  transporterNew(&T0, u.f, si0, &all, nelemV, nc);
  transporterNew(&Tb0, ua.f, sib0, &all, nelemV, nc);
  //transporterNew(&M0, NULL, si0, &all, nelemM, nc);
  transporterNew(&M0, u.f, si0, &all, nelemM, nc);

  int nrep = (int)(1 + (nrepfac*nthreads)/vol);
  //nrep = 20000;
  printf0("nrep: %i\n", nrep);
  t0 = dtime();
  for(int i=0; i<nrep; i++) {
    transporterDo(&T0, v2.f, v1.f, 2);
  }
  t0 = dtime() - t0;
  double flops = 8*nc*nc*vol;
  printf0("xportV: %10g %10g mf: %g\n",t0,1e6*t0/nrep,1e-6*flops*nrep/t0);
  t0 = dtime();
  for(int i=0; i<nrep; i++) {
    transporterDo(&M0, u2.f, u1.f, 2);
  }
  t0 = dtime() - t0;
  flops = 8*nc*nc*nc*vol;
  printf0("xportM: %10g %10g mf: %g\n",t0,1e6*t0/nrep,1e-6*flops*nrep/t0);
  beginParallelSubset(all) {
    F_eq_r(&v2, 0);
    F_eq_r(&v3, 0);
    F_eq_r(&u2, 0);
  } endParallelSubset;
  transporterDo(&T0, v2.f, v1.f, 2);
  transporterDo(&M0, u2.f, u1.f, 2);
  printf0("%-20.16g\t%-20.16g\n", nv1, nu1);
  beginParallelSubset(all) {
    double r[2];
    r[0] = lnorm2_F(&v2);
    r[1] = lnorm2_F(&u2);
    sumDoubleArray(r, 2);
    singleThread { nv1 = r[0]; nu1 = r[1]; }
  } endParallelSubset;
  printf0("%-20.16g\t%-20.16g\n", nv1, nu1);
  transporterDo(&Tb0, v3.f, v2.f, 2);
  beginParallelSubset(all) {
    X_vpeq_r_times_X(v3.f, -1, v1.f, v1.nelem, ti0, ti1);
    double r;
    r = lnorm2_F(&v3);
    sumDoubleArray(&r, 1);
    singleThread { nv1 = r; }
  } endParallelSubset;
  printf0("%-20.16g\n", nv1);

  transporterDo(&M0, u3.f, u1.f, 0);
  qopTransportM(&layout, u4.f, u.f, u1.f, dir, 1);
  beginParallelSubset(all) {
    //F_eq_r(&u4, 0);
    //P(X_peq_M_times_X)(u4.f, u.f, u3.f, nelemM, nc, 2, ti0, ti1);
    X_vpeq_r_times_X(u4.f, -1, u2.f, nelemM, ti0, ti1);
    double r;
    r = lnorm2_F(&u4);
    sumDoubleArray(&r, 1);
    singleThread { nv1 = r; }
  } endParallelSubset;
  printf0("%-20.16g\n", nv1);

  qopStagDslashFini();
  return 0;
}
