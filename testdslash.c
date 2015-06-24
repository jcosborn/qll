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
	vecr nb0 = NEG(b0);
	vecr na0 = NEG(a0);
	ST(M(ic,jc,0), nb0);
	ST(M(ic,jc,1), b1);
	ST(M(jc,ic,0), na0);
	ST(M(jc,ic,1), a1);
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
getplaq(real *u[], Layout *l, int nc)
{
  int nd = l->nDim;
  int nelemM = 2*nc*nc;
  Field uf[nd], *ufp[nd];
  for(int i=0; i<nd; i++) {
    fieldSet(&uf[i], u[i], l, nelemM);
    ufp[i] = &uf[i];
  }
  int np = (nd*(nd-1))/2;
  double p[np];
  gaugePlaq(p, ufp);
  for(int i=0; i<np; i++) {
    printf0("%i\t%g\n", i, p[i]);
  }
}

void
testDslash(real *v1, real *v2, real *tv, real *u[], real *u3[], double mass,
	   int nrep, char *sub, P(StaggeredLinks) *sl)
{
  int vol = layout.nSites;
  double flops = 72*8*vol;
#ifdef NAIK
  flops *= 2;
#endif
  if(sub[0]=='a') flops += 6*vol;
  else flops *= 0.5;
  printf0("testing Dslash sub: %s\n", sub);
  printf0("nrep: %i\n", nrep);
  // warmup
  for(int i=0; i<10; i++) {
    P(dslash2)(sl, v1, mass, v2);
  }
  QMP_barrier();
  P(dslashResetTimer)();
  double t0 = dtime();
  for(int i=0; i<nrep; i++) {
    P(dslash2)(sl, v1, mass, v2);
  }
  t0 = dtime() - t0;
  //t0 = t0/nrep;
  printf0("dslash: %10g %10g mf: %g\n", t0, 1e6*t0/nrep, 1e-6*flops*nrep/t0);
  P(dslashPrintTimer)(1.0/nrep);

  printf0("checking Dslash %s\n", sub);
  setV(v1);
  P(dslash2)(sl, v1, mass, v2);
  setV(tv);
  qopDslash(&layout, tv, u, u3, mass, v2, sub);
  compare(&layout, tv, v1, 6, 0);
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
  double nrepfac = 1e7;
  for(int i=0; i<nd; i++) ls[i] = 4;
  for(int i=1,j=0; argv[i]; i++) {
    if(argv[i][0]=='r') nrepfac = atof(argv[i]+1);
    else ls[j++] = atoi(argv[i]);
  }
  int nc = 3;
  int nelemV = 2*nc;
  int nelemM = 2*nc*nc;

  layout.nranks = QMP_get_number_of_nodes();
  layout.myrank = QMP_get_node_number();
  myrank = layout.myrank;
  printf0("qopStagDslashInit time: %g\n", t0);
  t0 = dtime();
  P(stagDslashSetup)(&layout, nd, ls, NULL);
  t0 = dtime() - t0;
  printf0("stagDslashSetup time: %g\n", t0);
  int vol = layout.nSites;
  int nthreads = 1;

  real *v1, *v2, *u[8], *uu3[8], *tm, *tv;
  real mass = 0.5;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }

  v1 = myalloc(nelemV*vol*sizeof(real));
  v2 = myalloc(nelemV*vol*sizeof(real));
  tv = myalloc(nelemV*vol*sizeof(real));
  for(int i=0; i<8; i++) {
    u[i] = myalloc(18*vol*sizeof(real));
    uu3[i] = myalloc(18*vol*sizeof(real));
  }
  tm = myalloc(18*vol*sizeof(real));
#pragma omp parallel for
  for(int i=0; i<nelemV*vol; i++) {
    v1[i] = 10000*layout.myrank + 100*(i/nelemV) + (i%nelemV);
    //v2[i] = 10000*layout.myrank + 100*(i/nelemV) + (i%nelemV);
    int t = 1 + layout.myrank + i;
    t = (t/64) + (t%64);
    t = (t/64) + (t%64);
    t = (t/64) + (t%64);
    t = (t/64) + (t%64);
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
      t = (t/64) + (t%64);
      t = (t/64) + (t%64);
      t = (t/64) + (t%64);
      t = (t/64) + (t%64);
      u[j][i] = t/64.;
      uu3[j][i] = (t - 32)/32.;
      //u[j][i] = (j==0) * (i==0);
      //uu3[j][i] = 0;
      //u[j][i] = j*18*vol + i;
      //uu3[j][i] = j*18*vol + i;
      //u[j][i] = 1000*(r+s) + 10*e + l;
      //uu3[j][i] = 1000*(r+s) + 10*e + l;
      //u[j][i] = cos(t);
      //uu3[j][i] = sin(t - 10);
      //u[j][i] = 2*((((i/layout.nSitesInner)%18)%8) == 0);
    }
  }
  norm2(v2, nelemV*vol);
  norm2(v2, 3*vol);
  getplaq(u, &layout, 3);

  //QMP_barrier();
  //TRACE_ALL;
  t0 = dtime();
  qopStagDslashSetup(&layout);
  t0 = dtime() - t0;
  printf0("qopStagDslashSetup time: %g\n", t0);
  //QMP_barrier();
  //TRACE_ALL;
  for(int i=0; i<4; i++) {
    int i2 = 2*i;
#if 1
    printf0("checking shift1[%i]\n", +(i+1)); fflush(stdout);
    //norm2(u[i2], 18*vol);
    P(shift)(u[i2+1], u[i2], 18, +(i+1), 1);
    qopShift(&layout, tm, u[i2], 18, +(i+1), 1);
    compare(&layout, tm, u[i2+1], 18, 0);
    printf0("checking shift3[%i]\n", +(i+1)); fflush(stdout);
    //norm2(uu3[i2], 18*vol);
    P(shift)(uu3[i2+1], uu3[i2], 18, +(i+1), 3);
    qopShift(&layout, tm, uu3[i2], 18, +(i+1), 3);
    compare(&layout, tm, uu3[i2+1], 18, 0);
#endif
    printf0("checking shift1[%i]\n", -(i+1)); fflush(stdout);
    P(shift)(u[i2+1], u[i2], 18, -(i+1), 1);
    qopShift(&layout, tm, u[i2], 18, -(i+1), 1);
    compare(&layout, tm, u[i2+1], 18, 0);
    printf0("checking shift3[%i]\n", -(i+1)); fflush(stdout);
    P(shift)(uu3[i2+1], uu3[i2], 18, -(i+1), 3);
    qopShift(&layout, tm, uu3[i2], 18, -(i+1), 3);
    compare(&layout, tm, uu3[i2+1], 18, 0);
    adjM(&layout, u[i2+1]);
    adjM(&layout, uu3[i2+1]);
  }
#ifdef NAIK
  real **u3 = uu3;
#else
  real **u3 = NULL;
#endif

  {
    Field fl[nd], ll[nd], *flp[nd], *llp[nd];
    Field fl2[nd], ll2[nd], *fl2p[nd], *ll2p[nd];
    Field uu[nd], ul[nd], *uup[nd], *ulp[nd];
    Subset all;
    layoutSubset(&all, &layout, "all");
    for(int i=0; i<nd; i++) {
      fieldNew(&fl[i], &layout, nelemM);
      fieldNew(&ll[i], &layout, nelemM);
      flp[i] = &fl[i];
      llp[i] = &ll[i];
      fieldNew(&fl2[i], &layout, nelemM);
      fieldNew(&ll2[i], &layout, nelemM);
      fl2p[i] = &fl2[i];
      ll2p[i] = &ll2[i];
      fieldSet(&uu[i], u[i], &layout, nelemM);
      fieldSet(&ul[i], u[i], &layout, nelemM);
      uup[i] = &uu[i];
      ulp[i] = &ul[i];
    }
    Fat7Coeffs coef = FAT7COEFFS_ZERO;
    coef.one_link = 0.3-((0.1*15)/16);
    coef.three_staple = 0.2/(4*6);   //  6
    coef.five_staple = 0.2/(16*24);  // 24
    coef.seven_staple = 0.2/(64*48); // 48
    coef.lepage = 0.1/(16*6);        //  6
    coef.naik = 0.2;
    smearFat7(flp, llp, &coef, uup, ulp);
    qopSmear(fl2p, ll2p, &coef, uup, ulp);
    for(int i=0; i<nd; i++) {
      printf0("comparing fl[%i]\n", i);
      compare(&layout, fl[i].f, fl2[i].f, nelemM, 1);
      printf0("comparing ll[%i]\n", i);
      compare(&layout, ll[i].f, ll2[i].f, nelemM, 1);
    }
    for(int i=0; i<nd; i++) {
      fieldFree(&fl[i]);
      fieldFree(&ll[i]);
      fieldFree(&fl2[i]);
      fieldFree(&ll2[i]);
    }
  }

  //qopStagDslashFini();
  //exit(1);

  P(StaggeredLinks) sla, sle, slo;
  P(createStaggeredLinks)(&sla, &layout, u, u3, "all");
  P(createStaggeredLinks)(&sle, &layout, u, u3, "even");
  P(createStaggeredLinks)(&slo, &layout, u, u3, "odd");

  P(StaggeredLinks) slfba, slfbe, slfbo;
  P(createStaggeredLinksFb)(&slfba, &layout, u, u3, "all");
  P(createStaggeredLinksFb)(&slfbe, &layout, u, u3, "even");
  P(createStaggeredLinksFb)(&slfbo, &layout, u, u3, "odd");

  int nrep = (int)(1 + (nrepfac*nthreads)/vol);

  testDslash(v1, v2, tv, u, u3, mass, nrep, "all", &sla);
  testDslash(v1, v2, tv, u, u3, mass, nrep, "all", &slfba);

  testDslash(v1, v2, tv, u, u3, mass, nrep, "even", &sle);
  testDslash(v1, v2, tv, u, u3, mass, nrep, "even", &slfbe);

  testDslash(v1, v2, tv, u, u3, mass, nrep, "odd", &slo);
  testDslash(v1, v2, tv, u, u3, mass, nrep, "odd", &slfbo);

  //char *sub = "all";
  char *sub = "even";
  //char *sub = "odd";

  P(StaggeredSolver) sse;
  P(createStaggeredSolver)(&sse, &layout, u, u3, "even");

  int its = 0;
  double rsq = 1e-6;
  int maxits = 1000;
  //nrep = 10;
  double flops0 = 72*8*vol;
#ifdef NAIK
  flops0 *= 2;
#endif
  flops0 += nelemV*vol;

  printf0("timing solve\n");
  nrep = 1 + (nrep/400);
  P(dslashResetTimer)();
  t0 = dtime();
  for(int i=0; i<nrep; i++) {
    //its = P(solve)(v1, u, u3, mass, v2, rsq, sub);
    its = P(solve2)(&sse, v1, mass, v2, rsq, maxits, sub);
  }
  t0 = dtime() - t0;
  double flops = flops0 + 36*vol;
  printf0("solve: %i %10g %10g mf: %g\n",its,t0,1e6*t0/nrep,its*1e-6*flops*nrep/t0);
  //dslashPrintTimer(1.0/nrep);

  mass = 0.5;
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
    xs[i] = myalloc(nelemV*vol*sizeof(real));
    ys[i] = myalloc(nelemV*vol*sizeof(real));
    setV(xs[i]);
    setV(ys[i]);
  }

  for(int im=1; im<=maxmasses; im++) {
    printf0("timing multishift solve masses: %i\n", im);
    P(dslashResetTimer)();
    t0 = dtime();
    for(int i=0; i<nrep; i++) {
      //its = P(solveMulti)(xs, u, u3, masses, im, v2, rsq, sub);
      its = P(solveMulti2)(&sse, xs, masses, im, v2, rsq, maxits, sub);
    }
    t0 = dtime() - t0;
    flops = flops0 + (24+12*im)*vol;
    printf0("solve %i: %i %10g %10g mf: %g\n",im,its,t0,1e6*t0/nrep,
	    its*1e-6*flops*nrep/t0);
    //dslashPrintTimer(1.0/nrep);
  }

  printf0("checking solve\n");
  setV(v1);
  norm2(v2, nelemV*vol);
  norm2(v2, 3*vol);
  //its = P(solve)(v1, u, u3, mass, v2, rsq, sub);
  its = P(solve2)(&sse, v1, mass, v2, rsq, maxits, sub);
  printf0("qll its: %i\n", its);
  setV(tv);
  //qopSolve(&layout, tv, u, u3, mass, v2, rsq, sub);
  qopSolve(&layout, tv, u, u3, mass, v2, rsq, "even");
  compare(&layout, tv, v1, nelemV, 1);

  int nmasses = 5;
  sub = "even";
  printf0("checking multishift solve\n");
  //its = P(solveMulti)(xs, u, u3, masses, nmasses, v2, rsq, sub);
  its = P(solveMulti2)(&sse, xs, masses, nmasses, v2, rsq, maxits, sub);
  printf0("qll its: %i\n", its);
  qopSolveMulti(&layout, ys, u, u3, masses, v2, nmasses, rsq, "even");
  for(int i=0; i<nmasses; i++) {
    printf0("comparing mass[%i] = %g\n", i, masses[i]);
    compare(&layout, ys[i], xs[i], nelemV, 1);
  }

  qopStagDslashFini();
  return 0;
}
