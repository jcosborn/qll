#include <common.h>

#define printnorm(x) {double n2=norm2_X(x,1,ti0,ti1); gsum(&n2,1,TC); if(TC->tid==0) {printf0("%s: %g\n", #x, n2);}}
#define printf00 if(TC->tid==0) printf

/* milti-shift CG */
static int
cgmsSolveT(P(Cgls) *cgls, real *out[], real *in, double shifts[], int nshifts,
	   P(Linop) *Aop, void *Aargs, P(Linop) *Mop, void *Margs,
	   double rsqin[], int itnlim, ThreadCtx *TC)
{
  double a[nshifts], b[nshifts], s[nshifts];
  double bo[nshifts], z[nshifts], zo[nshifts], zn[nshifts];
  double rsq, oldrsq, pkp;
  double insq;
  double rsqstop;
  int iteration=0;

  int nelem = cgls->nelem;
  int ti0 = ((TC->tid*(nelem))/TC->nid);
  int ti1 = (((TC->tid+1)*(nelem))/TC->nid);

  r_eq_norm2_V(&insq, in);
  gsum(&insq, 1, TC);
  if(insq==0) {
    for(int i=0; i<nshifts; i++) {
      V_eq_zero(out[i]);
    }
    return 0;
  }

  LinopInfo li = LINOP_INFO_DEFAULT;
  li.wantAnorm2 = 1;
#define A(a,b) Aop(a,b,&li,Aargs,TC)

  real *r = (real *) cgls->r1;
  //real *z = (real *) cgls->z;
  //real *p = (real *) cgls->p;
  real *Ap = (real *) cgls->Ap;
  real *ps[nshifts];
  ps[0] = cgls->p;
  for(int i=1; i<nshifts; i++) {
    ps[i] = myalloc(nelem*sizeof(vec));
  }

  V_eq_V(r, in);
  for(int i=0; i<nshifts; i++) {
    V_eq_zero(out[i]);
    V_eq_V(ps[i], r);
    zo[i] = z[i] = 1;
    bo[i] = -1;
    a[i] = 0;
    s[i] = 1;
  }

  rsqstop = rsqin[0] * insq;
  rsq = insq;
  //relnorm2 = 0;
  if(TC->tid==0 && cgls->verbose>0) {
    //printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", -1, rsq);
    //printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", iteration, rsq/insq);
    //VERB(LOW, "CGMS: rsqstop = %g\n", rsqstop);
  }

  while(1) {
    oldrsq = rsq;

    A(Ap, ps[0]);
    //if(shifts[imin]!=0.0) V_peq_r_times_V(Ap, shifts+imin, p, subset);
    iteration++;

    //r_eq_re_V_dot_V(&pkp, ps[0], Ap);
    //gsum(&pkp, 1, TC);
    pkp = li.Anorm2;
    pkp *= s[0]*s[0];
    if(pkp<=0) break;  // loss of precision in calculating pkp

    b[0] = rsq / pkp;
    zn[0] = 1;
    for(int i=1; i<nshifts; i++) {
      double c1;
      zn[i] = z[i]*zo[i]*bo[0];
      c1 = b[0]*a[0]*(zo[i]-z[i]);
      c1 += zo[i]*bo[0]*(1+shifts[i]*b[0]);
      if(c1!=0.0) zn[i] /= c1;
      else zn[i] = 0;
      if(z[i]!=0.0) b[i] = b[0]*zn[i]/z[i];
      else zn[i] = b[i] = 0;
    }

    for(int i=0; i<nshifts; i++) {
      double t = b[i] * s[i];
      //printf00("t: %g\n", t);
      //printnorm(ps[i]);
      V_peq_r_times_V(out[i], t, ps[i]);
      //printnorm(out[i]);
    }

    {
      double t = b[0] * s[0];
      V_meq_r_times_V(r, t, Ap);
    }
    r_eq_norm2_V(&rsq, r);
    gsum(&rsq, 1, TC);

#if 0
    /* compute FNAL norm if requested */
    /* here we look at the largest shift, since the FNAL norm is
       most stringent for that case */
    if(rsqin[imax]->relmin > 0) {
      t = b[imax] * s[imax];
      v_meq_r_times_v(r, t, Ap, n);
      relnorm2 = relnorm2_v(r, out[imax], n);
    }
#endif

    if(TC->tid==0 && cgls->verbose>0) {
      //VERB(HI, "CGMS: iter %i rsq = %g rel = %g\n", iteration, rsq,
      //   relnorm2);
      printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", iteration, rsq/insq);
    }

    if( (iteration>=itnlim) ||
        (rsqstop <= 0 || rsq<rsqstop) ) {
   //&&(res_arg[imax]->relmin <= 0 || relnorm2<res_arg[imax]->relmin)) ) {
      /* only way out */
      break;
    }

    a[0] = rsq / oldrsq;
    for(int i=1; i<nshifts; i++) {
      double c2 = z[i]*b[0];
      if(c2!=0.0) a[i] = a[0]*zn[i]*b[i]/c2;
      else a[i] = 0;
    }

    //s[0] *= a[0];
    {
      //double t = 1.0/s[0];
      //V_peq_r_times_V(ps[0], t, r);
      V_eq_scale_plus_V(ps[0], a[0], r);
    }
    for(int i=1; i<nshifts; i++) {
      //s[i] *= a[i];
      if(s[i]!=0) {
	//double t = zn[i]/s[i];
	//V_peq_r_times_V(ps[i], t, r);
	V_eq_scale_plus_r_times_V(ps[i], a[i], zn[i], r);
	//printf00("t: %g\n", t);
	//printf00("zn: %g\n", zn[i]);
	//printf00("s: %g\n", s[i]);
	//printnorm(r);
	//printnorm(ps[i]);
      }
    }

    for(int i=0; i<nshifts; i++) {
      bo[i] = b[i];
      zo[i] = z[i];
      z[i] = zn[i];
    }
  }

  for(int i=1; i<nshifts; i++) {
    free(ps[i]);
  }

#if 0
  for(i=0; i<nshifts; i++) {
    res_arg[i]->final_rsq = rsq/insq;
    res_arg[i]->final_rel = relnorm2;
    res_arg[i]->final_iter = iteration;
    res_arg[i]->final_restart = 0;
  }
#endif
  //VERB(LOW, "CGMS: done: iter %i rsq = %g\n", iteration, rsq);

  return iteration;
}

int
P(cgmsSolve)(P(Cgls) *cgls, real *x[], real *b, double shifts[], int nshifts,
	     P(Linop) *Aop, void *Aargs, P(Linop) *Mop, void *Margs,
	  double rsqs[], int itnlim)
{
  int its = 0;
  BEGIN_PARALLEL
    {
      int it = cgmsSolveT(cgls, x, b, shifts, nshifts, Aop, Aargs,
			  Mop, Margs, rsqs, itnlim, TC);
      if(TC->tid==0) its = it;
    }
  END_PARALLEL;
  return its;
}
