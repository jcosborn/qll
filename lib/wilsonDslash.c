#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include <omp.h>
#include "common.h"

//#define DeqRtimesD(x,r,y) X_veq_r_times_X(x,r,y,24,ti0,ti1)
#define DeqRtimesD(x,r,y) if((r)==0) { X_veq_r(x,0,24,ti0,ti1); }	\
 else { X_veq_r_times_X(x,r,y,24,ti0,ti1); }
#define DpeqRtimesD(x,r,y) X_vpeq_r_times_X(x,r,y,24,ti0,ti1)
#define DpeqSprojMtimesVXp(x,m,i,p,y,d,s) P(D_peq_spinproj_M_times_xDp)(x,m,i,p,y,d,s,ti0,ti1)
#define DmeqSprojMtimesVXp(x,m,i,p,y,d,s) P(D_meq_spinproj_M_times_xDp)(x,m,i,p,y,d,s,ti0,ti1)

static double dtstartrecv=0,dtsend=0,dtpack=0,dtbarsend=0,dtstartsend=0,dtmass=0,dtlocal=0,dtwaitsend=0,dtrecv=0,dtwaitrecv=0,dtbarrecv=0,dtblend=0,dtdonerecv=0;

#define STAMP(v) v = rdtsc()
#define DSTAMP(v) double v = rdtsc()
#define DSTAMPD(v,v0) double v = rdtsc() - v0
//#define STAMP(v) double v = 0; if0{ v = rdtsc(); }

static void
dslash2T(P(WilsonLinks) *d, real *x, real mass, real *y, ThreadCtx *TC)
{
  SUBSET_DEFS2;
  int tid = TC->tid;
#define SI(i) d->t[i].si
#define SB(i) d->t[i].sb
#define SBUF(i) SB(i).sbuf
#define RBUF(i) SB(i).rbuf
  //int nsi = d->n;

  DSTAMP(t0);
#if 0
  if(tid==0) {
    for(int i=0; i<nsi; i++) {
      if(SI(i).nRecvRanks>0) {
	startRecvBuf(&SB(i));
      }
    }
  }

  DSTAMP(t1);
  double tp=0,tbs=0,ts=0;
  for(int i=0; i<nsi; i++) {
    if(SI(i).nSendRanks>0) {
      DSTAMP(tt0);
      if(i==0) {
        setOff(SI(i).nSendSites);
      } else {
        setOffC(SI(i).nSendSites);
      }
      P(X_veq_pack_xX)(SBUF(i), SI(i).pack, SI(i).sendSites, y, 6, off, len);
      DSTAMP(tt1);
      T0WAIT;
      DSTAMP(tt2);
      if(tid==0) {
	startSendBuf(&SB(i));
	//printf("vvs: %i\tpack: %i\n", s->vvs, s->pack);
	//for(int i=0; i<8; i++) {
	//printf("idxSend[%i]: %i\n", i, s->idxSend[i]);
	//}
      }
      DSTAMP(tt3);
      tp += tt1-tt0;
      tbs += tt2-tt1;
      ts += tt3-tt2;
    }
  }
#endif

  DSTAMP(t2);
  setOffSubC2(d->t[0].sub);
  DeqRtimesD(x, mass, y);
  DSTAMP(t21);

  int sti0 = ti0;
  int sti1 = ti1;
  int block = 999999;
  for(int ti0=sti0; ti0<sti1; ti0+=block) {
    ti1 = ti0 + block;
    if(ti1>sti1) ti1 = sti1;
#if 1
  P(D_meq_spinproj0m_M_times_xDp)(x,d->t[0].u,SI(0).pidx,SI(0).perm,y,ti0,ti1);
  P(D_meq_spinproj0p_M_times_xDp)(x,d->t[1].u,SI(1).pidx,SI(1).perm,y,ti0,ti1);
  P(D_meq_spinproj1m_M_times_xDp)(x,d->t[2].u,SI(2).pidx,SI(2).perm,y,ti0,ti1);
  P(D_meq_spinproj1p_M_times_xDp)(x,d->t[3].u,SI(3).pidx,SI(3).perm,y,ti0,ti1);
  P(D_meq_spinproj2m_M_times_xDp)(x,d->t[4].u,SI(4).pidx,SI(4).perm,y,ti0,ti1);
  P(D_meq_spinproj2p_M_times_xDp)(x,d->t[5].u,SI(5).pidx,SI(5).perm,y,ti0,ti1);
  P(D_meq_spinproj3m_M_times_xDp)(x,d->t[6].u,SI(6).pidx,SI(6).perm,y,ti0,ti1);
  P(D_meq_spinproj3p_M_times_xDp)(x,d->t[7].u,SI(7).pidx,SI(7).perm,y,ti0,ti1);
#else
  int nsi = d->n;
  for(int i=0; i<nsi; i++) {
    int mu = i/2;
    int s = 1 - 2*(i&1);
    DmeqSprojMtimesVXp(x, d->t[i].u, SI(i).pidx, SI(i).perm, y, mu, -s);
  }
#endif
  }
  ti0 = sti0;
  ti1 = sti1;

  DSTAMP(t3);
  double tws=0,twr=0,tbr=0,tb=0,tdr=0;
#if 0
  for(int i=0; i<nsi; i++) {
    if(SI(i).nRecvDests>0) {
      DSTAMP(tt0);
      double tt1;
      if(SI(i).nRecvRanks>0) {
        if(tid==0) {
          waitRecvBuf(&SB(i));
        }
        STAMP(tt1);
        TWAIT0;
      } else {
        STAMP(tt1);
      }
      DSTAMP(tt2);
      if(SB(i).offr[tid]<0) {
	setOffSubC(d->t[i].sub);
	int i0=0; while(i0<SI(i).nRecvDests && SI(i).recvDests[i0]<ti0) i0++;
	int i1=i0; while(i1<SI(i).nRecvDests && SI(i).recvDests[i1]<ti1) i1++;
	SB(i).offr[tid] = i0;
	SB(i).lenr[tid] = i1;
        //printf0("%i %i: i0: %i\ti1: %i\trd0: %i\trd1: %i\n", i, tid, i0, i1, SI(i).recvDests[i0], (i1-1)<0?0:SI(i).recvDests[i1-1]);
      }
      ti0 = SB(i).offr[tid];
      ti1 = SB(i).lenr[tid];
      VXpeqMXtimesVXb(SI(i).recvDests, x, d->t[i].u, SI(i).blend,
		      SI(i).recvLocalSrcs, y, SI(i).recvRemoteSrcs, RBUF(i));
      DSTAMP(tt3);
      if(SI(i).nRecvRanks>0) {
	if(tid==0) {
	  doneRecvBuf(&SB(i));
	}
      }
      DSTAMP(tt4);
      twr += tt1-tt0;
      tbr += tt2-tt1;
      tb += tt3-tt2;
      tdr += tt4-tt3;
    }
    if(SI(i).nSendRanks>0) {
      DSTAMP(tt0);
      if(tid==0) {
	waitSendBuf(&SB(i));
      }
      DSTAMP(tt1);
      tws += tt1-tt0;
    }
  }
#endif

  DSTAMP(t4);
  if(tid==0) {
    double t1 = t0;
    double tp = 0;
    double tbs = 0;
    double ts = 0;
    dtstartrecv += t1 - t0;
    dtsend += t2 - t1;
    dtpack += tp;
    dtbarsend += tbs;
    dtstartsend += ts;
    dtmass += t21 - t2;
    dtlocal += t3 - t21;
    dtrecv += t4 - t3;
    dtwaitrecv += twr;
    dtbarrecv += tbr;
    dtblend += tb;
    dtdonerecv += tdr;
    dtwaitsend += tws;
  }
#undef SI
#undef SB
#undef SBUF
#undef RBUF
}

void
P(createWilsonLinks)(P(WilsonLinks) *sl, Layout *l, real *u[], char *sub)
{
  int nd = l->nDim;
  int nl = 2*nd;
  sl->n = nl;
  sl->t = myalloc(nl*sizeof(Transporter));
  sl->massFac = 1;
  sl->fb = 0;
  ShiftIndices *si = P(getStaggeredSI)(0, sub);
  if(sub[0]=='e' || sub[0]=='o') {
    sl->massFac = 0;
  }

  for(int i=0; i<nl; i++) {
    sl->t[i].u = u[i];
    sl->t[i].si = si[i];
    prepareShiftBuf(&(sl->t[i].sb), &si[i], 24*sizeof(real));
    layoutSubset(&(sl->t[i].sub), l, sub);
    sl->t[i].backward = 0;
  }
}

#if 0
void
P(createWilsonLinksFb)(P(WilsonLinks) *sl, Layout *l, real *u[], char *sub)
{
  int nd = l->nDim;
  int nl = 2*nd;
  if(u3) nl *= 2;
  sl->n = nl;
  sl->t = myalloc(nl*sizeof(Transporter));
  sl->massFac = 1;
  sl->fb = 1;
  ShiftIndices *si = P(getWilsonSI)((u3!=NULL), sub);
  if(sub[0]=='e' || sub[0]=='o') {
    sl->massFac = 0;
  }

  for(int i=0; i<nl; i++) {
    if(u3) {
      if((i&2)==0) { // forward
	if((i&1)==0) {
	  sl->t[i].u = u[i/2];
	} else {
	  sl->t[i].u = u3[i/2];
	}
	sl->t[i].backward = 0;
      } else {
	if((i&1)==0) {
	  sl->t[i].u = u[i/2-1];
	} else {
	  sl->t[i].u = u3[i/2-1];
	}
	sl->t[i].backward = 1;
      }
    } else {
      if((i&1)==0) { // forward
	sl->t[i].u = u[i];
	sl->t[i].backward = 0;
      } else {
	sl->t[i].u = u[i-1];
	sl->t[i].backward = 1;
      }
    }
    sl->t[i].si = si[i];
#ifdef SHIFT13
    if(u3) {
      if((i&1)==1) {
	ShiftBuf *sbp[2] = {&(sl->t[i-1].sb),&(sl->t[i].sb)};
	ShiftIndices *sip[2] = {&si[i-1], &si[i]};
	prepareShiftBufs(sbp, sip, 2, 6*sizeof(real));
      }
    } else
#endif
      prepareShiftBuf(&(sl->t[i].sb), &si[i], 6*sizeof(real));
    layoutSubset(&(sl->t[i].sub), l, sub);
  }
}
#endif

void
P(freeWilsonLinks)(P(WilsonLinks) *sl)
{
  int n = sl->n;
  for(int i=0; i<n; i++) {
    freeShiftBuf(&(sl->t[i].sb));
  }
  free(sl->t);
}

void
P(wilsonDslash)(P(WilsonLinks) *d, real *x, real mass, int sign, real *y)
{
  mass = d->massFac*(4+mass);
  BEGIN_PARALLEL
    {
      if(d->fb) {
	//dslashFb2T(d, x, mass, y, TC);
      } else {
	dslash2T(d, x, mass, y, TC);
      }
    }
  END_PARALLEL;
}

static void
dslashp2(real *Ax, real *x, LinopInfo *li, void *args, ThreadCtx *TC)
{
  SUBSET_DEFS2;
  P(WilsonSolver) *ss = (P(WilsonSolver)*)args;
  OMP_BARRIER;
  //if(tid==0) printf0("%p\t%i\t%i\n", x, da->sso->begin, da->sso->end);
  dslash2T(&ss->slo, ss->v, 0, x, TC);
  OMP_BARRIER;
  dslash2T(&ss->sle, Ax, 0, ss->v, TC);
  //OMP_BARRIER;
  setOffSubC2(ss->sube);
  DeqRtimesD(Ax, -1, Ax);
  double m = ss->mass;
  double m2 = m*m;
  DpeqRtimesD(Ax, m2, x);
  //setOffSub(*(da->sse));
  //VeqRtimesV(Ax, da->mass, x);
  if(li->wantAnorm2) {
    double nrm2[2];
    //setOffSub(*(da->sse));
    r_veq_norm2_X(nrm2, x, 6, ti0, ti1);
    setOffSubC2(ss->subo);
    r_veq_norm2_X(nrm2+1, ss->v, 6, ti0, ti1);
    gsum(nrm2, 2, TC);
    li->Anorm2 = m2 * nrm2[0] + nrm2[1];
  } else {
    OMP_BARRIER;
  }
}

void
P(createWilsonSolver)(P(WilsonSolver) *ss, Layout *l,
			 real *u[], char *sub)
{
  P(createWilsonLinks)(&ss->sle, l, u, "even");
  P(createWilsonLinks)(&ss->slo, l, u, "odd");
  int nelem = 24*l->nEvenOuter;
  ss->v = myalloc(24*l->nSites*sizeof(real));
  P(cglsInit)(&ss->cgls, nelem);
  layoutSubset(&ss->sube, l, "even");
  layoutSubset(&ss->subo, l, "odd");
  ss->noscale = 0;
}

void
P(freeWilsonSolver)(P(WilsonSolver) *ss)
{
  P(cglsFree)(&ss->cgls);
  free(ss->v);
  P(freeWilsonLinks)(&ss->slo);
  P(freeWilsonLinks)(&ss->sle);
}

int
P(wilsonSolve)(P(WilsonSolver) *ss, real *x, real mass, real *y,
	       double rsq, int maxits, char *sub)
{
  ss->mass = mass;
  if(rsq<0) {
    ss->cgls.verbose = 3;
    rsq = -rsq;
  }
  int its = P(cgSolve)(&ss->cgls, x, y, dslashp2, ss, NULL, NULL, rsq, maxits);
  ss->its = its;

  if(!ss->noscale) {
    BEGIN_PARALLEL
      {
	SUBSET_DEFS2;
	setOffSub2(ss->sube);
	DeqRtimesD(x, mass, x);
      }
    END_PARALLEL;
  }

  return its;
}

#if 0
int
P(solveMulti2)(P(WilsonSolver) *ss, real *x[], double masses[], int nmasses,
	       real *y, double rsq, int maxits, char *sub)
{
  ss->mass = masses[0];
  if(rsq<0) {
    ss->cgls.verbose = 3;
    rsq = -rsq;
  }
  double m2[nmasses];
  double rsqs[nmasses];
  for(int i=0; i<nmasses; i++) {
    m2[i] = masses[i]*masses[i] - masses[0]*masses[0];
    rsqs[i] = rsq;
  }
  int its = P(cgmsSolve)(&ss->cgls, x, y, m2, nmasses, dslashp2, ss,
			 NULL, NULL, rsqs, maxits);
  ss->its = its;

  if(!ss->noscale) {
    BEGIN_PARALLEL
      {
	SUBSET_DEFS;
	setOffSub(ss->sube);
	for(int i=0; i<nmasses; i++) {
	  VeqRtimesV(x[i], masses[i], x[i]);
	}
      }
    END_PARALLEL;
  }

  return its;
}
#endif

void
P(wilsonDslashResetTimer)(void)
{
#define DO(x) x = 0
  DO(dtstartrecv);
  DO(dtsend);
  DO(dtpack);
  DO(dtbarsend);
  DO(dtstartsend);
  DO(dtmass);
  DO(dtlocal);
  DO(dtrecv);
  DO(dtwaitsend);
  DO(dtwaitrecv);
  DO(dtbarrecv);
  DO(dtblend);
  DO(dtdonerecv);
#undef DO
}

void
P(wilsonDslashPrintTimer)(double scale)
{
  double dt = 0;
#define DO(x) dt += x
  DO(dtstartrecv);
  DO(dtsend);
  //DO(dtpack);
  //DO(dtbarsend);
  //DO(dtstartsend);
  DO(dtmass);
  DO(dtlocal);
  DO(dtrecv);
  //DO(dtwaitsend);
  //DO(dtwaitrecv);
  //DO(dtbarrecv);
  //DO(dtblend);
#undef DO
#define DO(x) printf0("%13s: %12.2f\t%6.2f\n", #x, scale*x, (100*x)/dt)
  DO(dtstartrecv);
  DO(dtsend);
  DO(dtpack);
  DO(dtbarsend);
  DO(dtstartsend);
  DO(dtmass);
  DO(dtlocal);
  DO(dtrecv);
  DO(dtwaitrecv);
  DO(dtbarrecv);
  DO(dtblend);
  DO(dtdonerecv);
  DO(dtwaitsend);
#undef DO
  printf0("%13s: %12.2f\n", "total", scale*dt);
}
