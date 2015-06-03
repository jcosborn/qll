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
compare(Layout *l, real *x, real *y, int nelem, int summary)
{
  int len = l->nSites*nelem;
  int errcount=0;
  //int imax = 0;
  double emax = 0;
  double sum2 = 0;
  for(int i=0; i<len; i++) {
    double d = x[i] - y[i];
    double e = fabs(d);
    if(e>0) {
      sum2 += e*e;
      if(e>emax) { emax = e; /*imax = i;*/ }
      if(!summary && errcount<10) {
	printf0("diff[%i] = %g\t%g\t%g\n", i, e, x[i], y[i]);
      }
      errcount++;
    }
  }
  QMP_sum_int(&errcount);
  QMP_sum_double(&sum2);
  double rms = sqrt(sum2/l->physVol);
  if(errcount>0) {
    printf0("found %i errors\n", errcount);
    printf0("max: %g\n", emax);
    printf0("rms: %g\n", rms);
  }
}

// negative adj
void
adjM(Layout *l, real *x)
{
  int nc = 3;
  int es = nc*nc*2*l->nSitesInner;
  for(int i=0; i<l->nSitesOuter; i++) {
    //vec (*m)[nc][2] = (vec (*)[nc][2]) &x[es*i];
    vec *m = (vec *) &x[es*i];
#define M(a,b,c) &m[(a*nc+b)*2+c]
    for(int ic=0; ic<nc; ic++) {
      for(int jc=ic; jc<nc; jc++) {
	vecr a0 = LD(M(ic,jc,0));
	vecr a1 = LD(M(ic,jc,1));
	vecr b0 = LD(M(jc,ic,0));
	vecr b1 = LD(M(jc,ic,1));
	//vecr nb0 = NEG(b0);
	//vecr na0 = NEG(a0);
	vecr nb1 = NEG(b1);
	vecr na1 = NEG(a1);
	ST(M(ic,jc,0), b0);
	ST(M(ic,jc,1), nb1);
	ST(M(jc,ic,0), a0);
	ST(M(jc,ic,1), na1);
      }
    }
  }
#undef M
}

void
setV(real *v)
{
#pragma omp parallel for
  for(int i=0; i<6*layout.nSites; i++) {
    v[i] = 10000*layout.myrank + 100*(i/6) + (i%6);
  }
}

void
norm2(real *x, int n)
{
  double s = 0;
  for(int i=0; i<n; i++) s += x[i]*x[i];
  printf0("norm2: %g\n", s);
}

void
testDslash(real *v1, real *v2, real *tv, real *u[], double mass,
	   int nrep, char *sub, P(WilsonLinks) *sl)
{
  int sign = 1;
  int vol = layout.nSites;
  double flops = (12+2*66+24)*8*vol;
  if(sub[0]=='a') flops += 24*vol;
  else flops *= 0.5;
  printf0("testing Dslash sub: %s\n", sub);
  printf0("nrep: %i\n", nrep);
  // warmup
  for(int i=0; i<10; i++) {
    P(wilsonDslash)(sl, v1, mass, sign, v2);
  }
  QMP_barrier();
  P(wilsonDslashResetTimer)();
  double t0 = dtime();
  for(int i=0; i<nrep; i++) {
    P(wilsonDslash)(sl, v1, mass, sign, v2);
  }
  t0 = dtime() - t0;
  //t0 = t0/nrep;
  printf0("dslash: %10g %10g mf: %g\n", t0, 1e6*t0/nrep, 1e-6*flops*nrep/t0);
  P(wilsonDslashPrintTimer)(1.0/nrep);

  printf0("checking Dslash %s\n", sub);
  setV(v1);
  P(wilsonDslash)(sl, v1, mass, sign, v2);
  setV(tv);
  qopWilsonDslash(&layout, tv, u, mass, sign, v2, sub);
  compare(&layout, tv, v1, 6, 0);
}

int
main(int argc, char *argv[])
{
  double t0 = dtime();
  qopWilsonDslashInit(&argc, &argv);
  t0 = dtime() - t0;

#define ND 4
  int nd = ND;
  int ls[ND];
  double nrepfac = 1e7;
  for(int i=0; i<nd; i++) ls[i] = 4;
  for(int i=1,j=0; argv[i]; i++) {
    if(argv[i][0]=='r') nrepfac = atof(argv[i]+1);
    else ls[j++] = atoi(argv[i]);
  }

  layout.nranks = QMP_get_number_of_nodes();
  layout.myrank = QMP_get_node_number();
  myrank = layout.myrank;
  printf0("qopStagDslashInit time: %g\n", t0);
  t0 = dtime();
  P(stagDslashSetup)(&layout, nd, ls, NULL);
  t0 = dtime() - t0;
  printf0("stagDslashSetup time: %g\n", t0);
  //qopWilsonDslashFini();
  //exit(0);
  int vol = layout.nSites;
  int nthreads = 1;

  real *v1, *v2, *u[8], *tm, *tv;
  real mass = 1;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }

  v1 = myalloc(24*vol*sizeof(real));
  v2 = myalloc(24*vol*sizeof(real));
  tv = myalloc(24*vol*sizeof(real));
  for(int i=0; i<8; i++) {
    u[i] = myalloc(18*vol*sizeof(real));
  }
  tm = myalloc(18*vol*sizeof(real));
#pragma omp parallel for
  for(int i=0; i<24*vol; i++) {
    v1[i] = 10000*layout.myrank + 100*(i/6) + (i%6);
    //v2[i] = 10000*layout.myrank + 100*(i/6) + (i%6);
    int t = 1 + layout.myrank + i;
    while(t>64) t = (t/64) + (t%64);
    v2[i] = t;
    //v2[i] = 1;
    //v2[i] = (layout.myrank+1)*(i==0);
  }
#pragma omp parallel for
  for(int i=0; i<18*vol; i++) {
    //int l = i%layout.nSitesInner;
    //int e = (i/layout.nSitesInner)%18;
    //int s = i/(18*layout.nSitesInner);
    int r = layout.myrank;
    for(int j=0; j<8; j++) {
      int t = 1 + j + r + i;
      while(t>64) t = (t/64) + (t%64);
      u[j][i] = t;
#if 0
      if(j==0) {
	u[j][i] = (e&1) * ((e/2)%3 == (e/6));
      } else {
	u[j][i] = 0*t;
      }
#endif
      //u[j][i] = (j==0) * (i==0);
      //u[j][i] = j*18*vol + i;
      //u[j][i] = cos(t);
    }
  }
  norm2(v2, 24*vol);
  norm2(v2, 12*vol);

  //QMP_barrier();
  //TRACE_ALL;
  t0 = dtime();
  qopWilsonDslashSetup(&layout);
  t0 = dtime() - t0;
  printf0("qopStagDslashSetup time: %g\n", t0);
  //QMP_barrier();
  //TRACE_ALL;
  for(int i=0; i<4; i++) {
    int i2 = 2*i;
    printf0("checking shift1[%i]\n", +(i+1)); fflush(stdout);
    //norm2(u[i2], 18*vol);
    P(shift)(u[i2+1], u[i2], 18, +(i+1), 1);
    qopShift(&layout, tm, u[i2], 18, +(i+1), 1);
    compare(&layout, tm, u[i2+1], 18, 0);
    printf0("checking shift1[%i]\n", -(i+1)); fflush(stdout);
    P(shift)(u[i2+1], u[i2], 18, -(i+1), 1);
    qopShift(&layout, tm, u[i2], 18, -(i+1), 1);
    compare(&layout, tm, u[i2+1], 18, 0);
    adjM(&layout, u[i2+1]);
  }

  for(int i=0; i<5; i++) {
    for(int j=0; j<2; j++) {
      printf0("checking proj[%i][%i]\n", i, j); fflush(stdout);
      int s = 1-2*(j&1);
      P(H_eq_spinproj_D)(v1, v2, i, s, 0, layout.nSitesOuter);
      qopSproj(&layout, tv, v2, i, s);
      compare(&layout, tv, v1, 12, 0);
      printf0("checking recon[%i][%i]\n", i, j); fflush(stdout);
      P(D_eq_spinrecon_H)(v1, v2, i, s, 0, layout.nSitesOuter);
      qopSrecon(&layout, tv, v2, i, s);
      compare(&layout, tv, v1, 12, 0);
    }
  }

  //exit(1);

  P(WilsonLinks) sla, sle, slo;
  P(createWilsonLinks)(&sla, &layout, u, "all");
  P(createWilsonLinks)(&sle, &layout, u, "even");
  P(createWilsonLinks)(&slo, &layout, u, "odd");

  //P(WilsonLinks) slfba, slfbe, slfbo;
  //P(createWilsonLinksFb)(&slfba, &layout, u, "all");
  //P(createWilsonLinksFb)(&slfbe, &layout, u, "even");
  //P(createWilsonLinksFb)(&slfbo, &layout, u, "odd");

  int nrep = (int)(1 + (nrepfac*nthreads)/vol);

  testDslash(v1, v2, tv, u, mass, nrep, "all", &sla);
  //testDslash(v1, v2, tv, u, u3, mass, nrep, "all", &slfba);

  testDslash(v1, v2, tv, u, mass, nrep, "even", &sle);
  //testDslash(v1, v2, tv, u, u3, mass, nrep, "even", &slfbe);

  testDslash(v1, v2, tv, u, mass, nrep, "odd", &slo);
  //testDslash(v1, v2, tv, u, u3, mass, nrep, "odd", &slfbo);

#if 0
  //char *sub = "all";
  char *sub = "even";
  //char *sub = "odd";

  P(WilsonSolver) sse;
  P(createWilsonSolver)(&sse, &layout, u, "even");

  double rsq = 1e-6;
  int maxits = 1000;
  //nrep = 10;
#ifdef NAIK
  double flops0 = (72*16+6)*vol;
#else
  double flops0 = (72*8+6)*vol;
#endif
  double flops = flops0;
  if(sub[0]!='a') flops *= 0.5;

  printf0("checking solve\n");
  mass = 10;
  setV(v1);
  norm2(v2, 24*vol);
  norm2(v2, 12*vol);
  int its = 0;
  its = P(wilsonSolve)(&sse, v1, mass, v2, rsq, maxits, sub);
  printf0("qll its: %i\n", its);
  setV(tv);
  qopWilsonSolve(&layout, tv, u, mass, v2, rsq, "even");
  compare(&layout, tv, v1, 24, 1);

  printf0("timing solve\n");
  nrep = 1 + (nrep/400);
  P(dslashResetTimer)();
  t0 = dtime();
  for(int i=0; i<nrep; i++) {
    its = P(wilsonSolve)(&sse, v1, mass, v2, rsq, maxits, sub);
  }
  t0 = dtime() - t0;
  flops = flops0 + 36*vol;
  printf0("solve: %10g %10g mf: %g\n", t0,1e6*t0/nrep,its*1e-6*flops*nrep/t0);
  //dslashPrintTimer(1.0/nrep);
#endif

#if 0
  int maxmasses = 10;
  double masses[maxmasses];
  masses[0] = mass;
  real *xs[maxmasses], *ys[maxmasses];
  xs[0] = v1;
  ys[0] = tv;
  setV(xs[0]);
  setV(ys[0]);
  for(int i=1; i<maxmasses; i++) {
    masses[i] = 1.4*masses[i-1];
    xs[i] = myalloc(6*vol*sizeof(real));
    ys[i] = myalloc(6*vol*sizeof(real));
    setV(xs[i]);
    setV(ys[i]);
  }

  for(int im=1; im<=maxmasses; im++) {
    printf0("timing multishift solve masses: %i\n", im);
    P(dslashResetTimer)();
    t0 = dtime();
    for(int i=0; i<nrep; i++) {
      //its = P(solveMulti)(xs, u, u3, masses, im, v2, rsq, sub);
      its = P(wilsonSolveMulti)(&sse, xs, masses, im, v2, rsq, maxits, sub);
    }
    t0 = dtime() - t0;
    flops = flops0 + (24+12*im)*vol;
    printf0("solve %i: %10g %10g mf: %g\n",im,t0,1e6*t0/nrep,
	    its*1e-6*flops*nrep/t0);
    //dslashPrintTimer(1.0/nrep);
  }

  int nmasses = 5;
  sub = "even";
  printf0("checking multishift solve\n");
  its = P(wilsonSolveMulti)(&sse, xs, masses, nmasses, v2, rsq, maxits, sub);
  printf0("qll its: %i\n", its);
  qopSolveMulti(&layout, ys, u, u3, masses, v2, nmasses, rsq, "even");
  for(int i=0; i<nmasses; i++) {
    printf0("comparing mass[%i] = %g\n", i, masses[i]);
    compare(&layout, ys[i], xs[i], 6, 1);
  }
#endif

  qopWilsonDslashFini();
  return 0;
}
