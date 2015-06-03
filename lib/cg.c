#include <stdlib.h>
#include <string.h>
#include "common.h"

#define create(v) cgls->v = myalloc(nelem*sizeof(vec))
#define destroy(v) free(cgls->v); cgls->v=NULL

void
P(cglsInit)(P(Cgls) *cgls, int nelem)
{
  cgls->nelem = nelem;
  cgls->verbose = 0;
  cgls->indent = 0;
  create(r1);
  create(r2);
  create(z);
  create(p);
  create(Ap);
  create(Adr2);
  //printf0("ap: %p\n", cgls->p);
}

void
P(cglsFree)(P(Cgls) *cgls)
{
  //printf0("dp: %p\n", cgls->p);
  destroy(r1);
  destroy(r2);
  destroy(z);
  destroy(p);
  destroy(Ap);
  destroy(Adr2);
}

#define seti(t) if(!strcmp(s,#t)) cgls->t = (int) v;
void
P(cglsSet)(P(Cgls) *cgls, char *s, double v)
{
  seti(verbose);
  seti(indent);
}

// solves  A x = b  for HPD A
static int
cgSolveT(P(Cgls) *cgls, real *x, real *b, P(Linop) *Aop, void *Aargs,
	    P(Linop) *Mop, void *Margs, double rsqin, int itnlim, ThreadCtx *TC)
{
  double bsq, rsq, rsqstop, alpha, beta, pAp;
  int itn = 0;
  int nelem = cgls->nelem;
  int ti0 = ((TC->tid*(nelem))/TC->nid);
  int ti1 = (((TC->tid+1)*(nelem))/TC->nid);
  LinopInfo li = LINOP_INFO_DEFAULT;
  li.wantAnorm2 = 1;

  real *r = (real *) cgls->r1;
  real *z = (real *) cgls->z;
  real *p = (real *) cgls->p;
  real *Ap = (real *) cgls->Ap;
  //if(TC->tid==0) printf0("p: %p\n", p);

#define A(a,b) Aop(a,b,&li,Aargs,TC)

  V_eq_zero(x);
  V_eq_zero(p);
  V_eq_V(r, b);

  r_eq_norm2_V(&rsq, r);
  gsum(&rsq, 1, TC);
  bsq = rsq;
  rsqstop = rsqin*bsq;

  if(TC->tid==0 && cgls->verbose>0) {
    printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", -1, rsq);
    printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", itn, rsq/bsq);
  }
  do {
    itn++;
    if(Mop) {
      Mop(z, r, &li, Margs, TC);
      r_eq_re_V_dot_V(&beta, r, z);
      gsum(&beta, 1, TC);
      beta = 1/beta;
      V_peq_r_times_V(p, beta, z);
    } else {
      beta = 1/rsq;
      V_peq_r_times_V(p, beta, r);
    }
    //r_eq_re_V_dot_V(&beta, p, r1);
    //printf("%g\n", beta);
    A(Ap, p);
    //r_eq_norm2_V(&Ap2, Ap);
    //r_eq_re_V_dot_V(&pAp, p, Ap);
    //gsum(&pAp, 1, TC);
    //if(TC->tid==0) printf0("pAp: %g  Anorm2: %g\n", pAp, li.Anorm2);
    pAp = li.Anorm2;
    alpha = 1/pAp;
    //QOP_printf0("alpha = %g\n", alpha);
    V_peq_r_times_V(x, alpha, p);
    V_meq_r_times_V(r, alpha, Ap);
#if 0
    QLA_Complex zz;
    c_eq_V_dot_V(&zz, p, r1);
    QOP_printf0("%g %g\n", QLA_real(zz), QLA_imag(zz));
    c_eq_V_dot_V(&zz, p, b1);
    QOP_printf0("%g %g\n", QLA_real(zz), QLA_imag(zz));
#endif

    r_eq_norm2_V(&rsq, r);
    gsum(&rsq, 1, TC);

    if(TC->tid==0 && cgls->verbose>1) {
      printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", itn, rsq/bsq);
    }
  } while(rsq>=rsqstop && itn<itnlim);
  if(TC->tid==0 && cgls->verbose>0) {
    printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", itn, rsq/bsq);
  }

  return itn;
#undef A
}

int
P(cgSolve)(P(Cgls) *cgls, real *x, real *b, P(Linop) *Aop, void *Aargs,
	   P(Linop) *Mop, void *Margs, double rsq, int itnlim)
{
  int its = 0;
  BEGIN_PARALLEL
    {
      int it = cgSolveT(cgls, x, b, Aop, Aargs, Mop, Margs, rsq, itnlim, TC);
      if(TC->tid==0) its = it;
    }
  END_PARALLEL;
  return its;
}
