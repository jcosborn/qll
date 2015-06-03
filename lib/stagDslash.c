#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include <omp.h>
#include "common.h"

//#define VeqRtimesV(x,r,y) P(V_veq_r_times_V)(x,r,y,off,len)
//#define VeqRtimesV(x,r,y) X_veq_r_times_X(x,r,y,6,ti0,ti1)
#define VeqRtimesV(x,r,y) if((r)==0) { X_veq_r(x,0,6,ti0,ti1); }	\
 else { X_veq_r_times_X(x,r,y,6,ti0,ti1); }

#define VpeqRtimesV(x,r,y) X_veq_X_plus_r_times_X(x,x,r,y,6,ti0,ti1)
#define VpeqMtimesVX(x,m,i,y) P(V_vpeq_M_times_xV)(x,m,i,y,off,len)
#define VpeqMtimesVXc(x,m,i,y) P(V_vpeq_M_times_xVc)(x,m,i,y,off,len)
#define VpeqMtimesVXp(x,m,i,p,y) P(V_vpeq_M_times_xVp)(x,m,i,p,y,ti0,ti1)
#define VmeqMAXtimesVXp(x,m,i,p,y) P(V_vmeq_xMap_times_xVp)(x,m,i,p,y,ti0,ti1)
#define VpeqMtimesVXn(x,m,i,y) P(V_vpeq_M_times_xVn)(x,m,i,y,off,len)
#define VpeqMtimesVXb(x,m,b,i1,y1,i2,y2) P(V_vpeq_M_times_xVb)(x,m,b,i1,y1,i2,y2,ti0,ti1)
#define VXpeqMXtimesVXb(i0,x,m,b,i1,y1,i2,y2) P(xV_vpeq_xM_times_yVzVb)(i0,x,m,b,i1,y1,i2,y2,ti0,ti1)
#define VXpeqMXtimesVY(ix,x,m,iy,y) P(xV_vpeq_xM_times_yV)(ix,x,m,iy,y,off,len)
#define VXpeqMYtimesVY(ix,x,iy,m,y) P(xV_vpeq_yM_times_yV)(ix,x,iy,m,y,off,len)
#define VXpeqMtimesV(ix,x,m,y) P(xV_vpeq_M_times_V)(ix,x,m,y,off,len)
#define VXpeqMXtimesV(ix,x,m,y) P(xV_vpeq_xM_times_V)(ix,x,m,y,off,len)

static void
shiftT(real *x, real *y, int nelem, ShiftIndices *s, ShiftBuf *sb,
       ThreadCtx *TC)
{
  SUBSET_DEFS2;
  int tid = TC->tid;

  if(s->nRecvRanks>0) {
    if(tid==0) {
      startRecvBuf(sb);
    }
  }
  if(s->nSendRanks>0) {
    //if(tid==0) {
    //printf0("nss: %i\tnss1: %i\tpack: %i\n", s->nSendSites, s->nSendSites1, s->pack);
    //for(int i=0; i<8; i++) {
    //printf0("sendSites[%i]: %i\n", i, s->sendSites[i]);
    //}
    //fflush(stdout);
    //QMP_barrier();
    //}
    //OMP_BARRIER;
    setOff2(s->nSendSites);
    P(X_eq_pack_xX)(sb->sbuf, s->pack, s->sendSites, y, nelem, ti0, ti1);
    OMP_BARRIER;
    if(tid==0) {
      startSendBuf(sb);
      //printf("vvs: %i\tpack: %i\n", s->vvs, s->pack);
      //for(int i=0; i<8; i++) {
      //printf("idxSend[%i]: %i\n", i, s->idxSend[i]);
      //}
    }
  }

  setOff2(s->vv);
  //X_veq_xXp(x, s->sidx, s->perm, s->idxSend, y, nelem, off, len);
  P(X_eq_xXp)(x, s->pidx, s->perm, y, nelem, ti0, ti1);

  if(s->nSendRanks>0) {
    if(tid==0) {
      waitSendBuf(sb);
    }
  }
  if(s->nRecvRanks>0) {
    if(tid==0) {
      //printf("vv: %i\tblend: %i\n", s->vv, s->blend);
      //for(int i=0; i<16; i++) {
      //printf("sidx[%i]: %i\n", i, s->sidx[i]);
      //}
      waitRecvBuf(sb);
    }
    OMP_BARRIER;
    //X_veq_blend_xXxXn(x, s->blend, s->sendSites,y,s->sidx,rbuf,nelem,off,len);
    setOff2(s->nRecvDests);
    int off = ti0;
    int len = ti1-ti0;
    P(xX_veq_blend_xXxX)(s->recvDests, x, s->blend, s->recvLocalSrcs, y,
			 s->recvRemoteSrcs, sb->rbuf, nelem, off, len);
  }
}

void
P(shift)(real *x, real *y, int nelem, int dir, int len)
{
  //int sd = 2*(abs(dir)-1);
  //if(dir<0) sd++;
  //ShiftIndices *s = (len==1) ? &si1[sd] : &si3[sd];
  ShiftIndices *s = P(getShift)(abs(dir)-1, (dir>0)?len:-len, "all");
  ShiftBuf sbs;
  prepareShiftBuf(&sbs, s, nelem*sizeof(real));
  BEGIN_PARALLEL
    {
      shiftT(x, y, nelem, s, &sbs, TC);
    }
  END_PARALLEL
  freeShiftBuf(&sbs);
}

static double dtstartrecv=0,dtsend=0,dtpack=0,dtbarsend=0,dtstartsend=0,dtmass=0,dtlocal=0,dtwaitsend=0,dtrecv=0,dtwaitrecv=0,dtbarrecv=0,dtblend=0,dtdonerecv=0;

#define STAMP(v) v = rdtsc()
#define DSTAMP(v) double v = rdtsc()
#define DSTAMPD(v,v0) double v = rdtsc() - v0
//#define STAMP(v) double v = 0; if0{ v = rdtsc(); }

static void
dslashFb2T(P(StaggeredLinks) *sl, real *x, real mass, real *y, ThreadCtx *TC)
{
  SUBSET_DEFS2;
  int tid = TC->tid;
  int nsi = sl->n;
  //int nn = nsi/2;
#define SI(i) sl->t[i].si
#define SB(i) sl->t[i].sb
#define SBUF(i) SB(i).sbuf
#define RBUF(i) SB(i).rbuf

  DSTAMP(t0);
#if 0
  if(tid==0) {
    for(int i=0; i<NSI; i++) {
      if(SI(i).nRecvRanks>0) {
	startRecvBuf(&SB(i));
      }
    }
  }

  DSTAMP(t1);
  double tp=0,tbs=0,ts=0;
  for(int i=0; i<NSI; i++) {
    if(SI(i).nSendRanks>0) {
      DSTAMP(tt0);
      if(i==0) {
        setOff2(SI(i).nSendSites);
      } else {
        setOffC2(SI(i).nSendSites);
      }
      //X_veq_pack_xX(SBUF(i), SI(i).pack, SI(i).sendSites, y, 6, off, len);
      X_eq_pack_xX(SBUF(i), SI(i).pack, SI(i).sendSites, y, 6, ti0, ti1);
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
  setOffSubC2(sl->t[0].sub);
  VeqRtimesV(x, mass, y);
  DSTAMP(t21);

  for(int i=0; i<nsi; i++) {
    if(sl->t[i].backward) {
      VmeqMAXtimesVXp(x, sl->t[i].u, SI(i).pidx, SI(i).perm, y);
    } else {
      VpeqMtimesVXp(x, sl->t[i].u, SI(i).pidx, SI(i).perm, y);
    }
  }

  DSTAMP(t3);

#if 0
  double tws=0,twr=0,tbr=0,tb=0,tdr=0;
  for(int i=0; i<NSI; i++) {
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
	setOffSubC(*sub);
	int i0=0; while(i0<SI(i).nRecvDests && SI(i).recvDests[i0]<ti0) i0++;
	int i1=i0; while(i1<SI(i).nRecvDests && SI(i).recvDests[i1]<ti1) i1++;
	SB(i).offr[tid] = i0;
	SB(i).lenr[tid] = i1;
        //printf0("%i %i: i0: %i\ti1: %i\trd0: %i\trd1: %i\n", i, tid, i0, i1, SI(i).recvDests[i0], (i1-1)<0?0:SI(i).recvDests[i1-1]);
      }
      ti0 = SB(i).offr[tid];
      ti1 = SB(i).lenr[tid];
      VXpeqMXtimesVXb(SI(i).recvDests,x,u[i],SI(i).blend,SI(i).recvLocalSrcs,
		      y,SI(i).recvRemoteSrcs,RBUF(i));
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

  double t1 = t0;
  double tp = 0;
  double tbs = 0;
  double ts = 0;
  double twr = 0;
  double tbr = 0;
  double tb = 0;
  double tdr = 0;
  double tws = 0;
  DSTAMP(t4);
  if(tid==0) {
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

static void
dslash2T(P(StaggeredLinks) *d, real *x, real mass, real *y, ThreadCtx *TC)
{
  SUBSET_DEFS;
  int tid = TC->tid;
#define SI(i) d->t[i].si
#define SB(i) d->t[i].sb
#define SBUF(i) SB(i).sbuf
#define RBUF(i) SB(i).rbuf
  int nsi = d->n;

  DSTAMP(t0);
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

  DSTAMP(t2);
  setOffSubC(d->t[0].sub);
  //printf("ti0: %i  ti1: %i\n", ti0, ti1); fflush(stdout);
  //printf("off: %i  len: %i\n", off, len); fflush(stdout);
  //printf("x: %p  y: %p\n", x, y);
  //printf("x[0]: %g  y[0]: %g\n", x[0], y[0]); fflush(stdout);
  VeqRtimesV(x, mass, y);
  DSTAMP(t21);

#if 1
  int sti0 = ti0;
  int sti1 = ti1;
  //int block = 32;
  int block = 256;
  //int block = 999999;
  for(int ti0=sti0; ti0<sti1; ti0+=block) {
    ti1 = ti0 + block;
    if(ti1>sti1) ti1 = sti1;
    for(int i=0; i<nsi; i++) {
      VpeqMtimesVXp(x, d->t[i].u, SI(i).pidx, SI(i).perm, y);
    }
  }
  ti0 = sti0;
  ti1 = sti1;
#else
  int sti0 = ti0;
  int sti1 = ti1;
  //int block = 1024*1024;
  int block = 8;
  int stride = block*TC->nid;
  int i0 = d->t[0].sub.beginOuter + block*TC->tid;
  int i1 = d->t[0].sub.endOuter;
  for(int ti0=i0; ti0<i1; ti0+=stride) {
    ti1 = ti0 + block;
    if(ti1>i1) ti1 = i1;
    for(int i=0; i<nsi; i++) {
      VpeqMtimesVXp(x, d->t[i].u, SI(i).pidx, SI(i).perm, y);
    }
  }
  ti0 = sti0;
  ti1 = sti1;
#endif

  DSTAMP(t3);
  double tws=0,twr=0,tbr=0,tb=0,tdr=0;
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

  DSTAMP(t4);
  if(tid==0) {
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
P(createStaggeredLinks)(P(StaggeredLinks) *sl, Layout *l,
			real *u[], real *u3[], char *sub)
{
  int nd = l->nDim;
  int nl = 2*nd;
  if(u3) nl *= 2;
  sl->n = nl;
  sl->t = myalloc(nl*sizeof(Transporter));
  sl->massFac = 1;
  sl->fb = 0;
  ShiftIndices *si = P(getStaggeredSI)((u3!=NULL), sub);
  if(sub[0]=='e' || sub[0]=='o') {
    sl->massFac = 0;
  }

  for(int i=0; i<nl; i++) {
    if(u3) {
      if((i&1)==0) {
	sl->t[i].u = u[i/2];
      } else {
	sl->t[i].u = u3[i/2];
      }
    } else {
      sl->t[i].u = u[i];
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
    sl->t[i].backward = 0;
  }
}

void
P(createStaggeredLinksFb)(P(StaggeredLinks) *sl, Layout *l,
			  real *u[], real *u3[], char *sub)
{
  int nd = l->nDim;
  int nl = 2*nd;
  if(u3) nl *= 2;
  sl->n = nl;
  sl->t = myalloc(nl*sizeof(Transporter));
  sl->massFac = 1;
  sl->fb = 1;
  ShiftIndices *si = P(getStaggeredSI)((u3!=NULL), sub);
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

void
P(freeStaggeredLinks)(P(StaggeredLinks) *sl)
{
  int n = sl->n;
  for(int i=0; i<n; i++) {
    freeShiftBuf(&(sl->t[i].sb));
  }
  free(sl->t);
}

void
P(dslash2)(P(StaggeredLinks) *d, real *x, real mass, real *y)
{
  mass *= d->massFac;
  BEGIN_PARALLEL
    {
      if(d->fb) {
	dslashFb2T(d, x, mass, y, TC);
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
  P(StaggeredSolver) *ss = (P(StaggeredSolver)*)args;
  OMP_BARRIER;
  //if(tid==0) printf0("%p\t%i\t%i\n", x, da->sso->begin, da->sso->end);
  dslash2T(&ss->slo, ss->v, 0, x, TC);
  OMP_BARRIER;
  dslash2T(&ss->sle, Ax, 0, ss->v, TC);
  //OMP_BARRIER;
  setOffSubC2(ss->sube);
  VeqRtimesV(Ax, -1, Ax);
  double m = ss->mass;
  double m2 = m*m;
  VpeqRtimesV(Ax, m2, x);
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
P(createStaggeredSolver)(P(StaggeredSolver) *ss, Layout *l,
			 real *u[], real *u3[], char *sub)
{
  P(createStaggeredLinks)(&ss->sle, l, u, u3, "even");
  P(createStaggeredLinks)(&ss->slo, l, u, u3, "odd");
  int nelem = 6*l->nEvenOuter;
  ss->v = myalloc(6*l->nSites*sizeof(real));
  P(cglsInit)(&ss->cgls, nelem);
  layoutSubset(&ss->sube, l, "even");
  layoutSubset(&ss->subo, l, "odd");
  ss->noscale = 0;
}

void
P(freeStaggeredSolver)(P(StaggeredSolver) *ss)
{
  P(cglsFree)(&ss->cgls);
  free(ss->v);
  P(freeStaggeredLinks)(&ss->slo);
  P(freeStaggeredLinks)(&ss->sle);
}

int
P(solve2)(P(StaggeredSolver) *ss, real *x, real mass, real *y,
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
	VeqRtimesV(x, mass, x);
      }
    END_PARALLEL;
  }

  return its;
}

int
P(solveMulti2)(P(StaggeredSolver) *ss, real *x[], double masses[], int nmasses,
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
	SUBSET_DEFS2;
	setOffSub2(ss->sube);
	for(int i=0; i<nmasses; i++) {
	  VeqRtimesV(x[i], masses[i], x[i]);
	}
      }
    END_PARALLEL;
  }

  return its;
}

void
P(dslashResetTimer)(void)
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
P(dslashPrintTimer)(double scale)
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

#define ND 4
#define ND2 (2*ND)
static ShiftIndices si1[ND2], si3[ND2], si13[2*ND2];
static ShiftIndices si1e[ND2], si3e[ND2], si13e[2*ND2];
static ShiftIndices si1o[ND2], si3o[ND2], si13o[2*ND2];

ShiftIndices *
P(getStaggeredSI)(int naik, char *sub)
{
  ShiftIndices *si;
  if(naik) si = si13; else si = si1;
  if(sub[0]=='e') {
    if(naik) si = si13e; else si = si1e;
  }
  if(sub[0]=='o') {
    if(naik) si = si13o; else si = si1o;
  }
  return si;
}

ShiftIndices *
P(getShift)(int dir, int len, char *sub)
{
  ShiftIndices *si;
  if(abs(len)==1) si = si1; else si = si3;
  if(sub[0]=='e') {
    if(abs(len)==1) si = si1e; else si = si3e;
  }
  if(sub[0]=='o') {
    if(abs(len)==1) si = si1o; else si = si3o;
  }
  if(len>0) si += 2*dir;
  else si += 2*dir + 1;
  return si;
}

// partition x into n blocks with geometry nx
// dist=0: prefer to split already split direction
// dist=1: split unsplit directions first
static void
partitionGeom(int *lx, int *nx, int *x, int n, int nd, int dist)
{
  for(int i=0; i<nd; i++) {
    nx[i] = 1;
    lx[i] = x[i];
  }
  int ww = n;
  while(ww>1) {
    int k = nd-1;
    while(lx[k]&1) {
      k--;
      if(k<0) {
	printf("not enough 2's in partitioned geom:");
	for(int i=0; i<nd; i++) printf(" %i", x[i]);
	printf(" /%i\n", n);
	exit(-1);
      }
    }
    for(int i=k-1; i>=0; i--) {
      if((lx[i]&1)==0) {
	if(dist==0) {
	  if(lx[i]>lx[k] || (lx[i]==lx[k] && nx[i]>nx[k])) k = i;
	} else {
	  if(nx[i]<nx[k] || (nx[i]==nx[k] && lx[i]>lx[k])) k = i;
	}
      }
    }
    nx[k] *= 2;
    lx[k] /= 2;
    ww /= 2;
  }
}

void
P(stagDslashSetup)(Layout *layout, int nd, int ls[], int rg0[])
{
  int rg[nd], lg[nd], og[nd], ig[nd];
  myrank = layout->myrank;
  PRINTV0("#physGeom:", "%i", ls, nd);
  if(rg0) {
    for(int i=0; i<nd; i++) {
      rg[i] = rg0[i];
      lg[i] = ls[i]/rg[i];
    }
  } else {
    partitionGeom(lg, rg, ls, layout->nranks, nd, 1);
  }
  PRINTV0("#rankGeom:", "%i", rg, nd);
  PRINTV0("#localGeom:", "%i", lg, nd);
  partitionGeom(og, ig, lg, VLEN, nd, 1);
  PRINTV0("#outerGeom:", "%i", og, nd);
  PRINTV0("#innerGeom:", "%i", ig, nd);

  layout->nDim = nd;
  layout->physGeom = malloc(nd*sizeof(int));
  layout->rankGeom = malloc(nd*sizeof(int));
  layout->innerGeom = malloc(nd*sizeof(int));
  for(int i=0; i<nd; i++) {
    layout->physGeom[i] = ls[i];
    layout->rankGeom[i] = rg[i];
    layout->innerGeom[i] = ig[i];
  }
  layoutSetup(layout);

  //TRACE_ALL;
  for(int i=0; i<8; i++) {
    int disp[2][nd];
    for(int j=0; j<nd; j++) disp[0][j] = 0;
    disp[0][i/2] = 2*(i&1) - 1;  // -ve is QDP_forward (from forward)
    for(int j=0; j<nd; j++) disp[1][j] = 0;
    disp[1][i/2] = 3*(2*(i&1) - 1);  // -ve is QDP_forward (from forward)

    makeShift(&si1[i], layout, disp[0]);
    makeShift(&si3[i], layout, disp[1]);

    makeShiftSub(&si1e[i], layout, disp[0], "even");
    makeShiftSub(&si3e[i], layout, disp[1], "even");

    makeShiftSub(&si1o[i], layout, disp[0], "odd");
    makeShiftSub(&si3o[i], layout, disp[1], "odd");

#ifdef SHIFT13
#if 0
    if((i&1)==0) {
      int dispb[2][nd];
      for(int j=0; j<nd; j++) dispb[0][j] = 0;
      dispb[0][i/2] = 2*(1) - 1;  // -ve is QDP_forward (from forward)
      for(int j=0; j<nd; j++) dispb[1][j] = 0;
      dispb[1][i/2] = 3*(2*(1) - 1);  // -ve is QDP_forward (from forward)

      ShiftIndices *si13p[4] = {&si13[2*i],&si13[2*i+1],&si13[2*(i+1)],&si13[2*(i+1)+1]};
      int *disps[4] = {disp[0],disp[1],dispb[0],dispb[1]};
      makeShiftMulti(si13p, layout, disps, 4);
    }
#else
    ShiftIndices *si13p[2] = {&si13[2*i],&si13[2*i+1]};
    int *disps[2] = {disp[0],disp[1]};
    makeShiftMulti(si13p, layout, disps, 2);

    si13p[0] = &si13e[2*i];
    si13p[1] = &si13e[2*i+1];
    char *subs[2] = {"even","even"};
    makeShiftMultiSub(si13p, layout, disps, subs, 2);

    si13p[0] = &si13o[2*i];
    si13p[1] = &si13o[2*i+1];
    subs[0] = "odd";
    subs[1] = "odd";
    makeShiftMultiSub(si13p, layout, disps, subs, 2);
#endif
#else // not SHIFT13
    si13[2*i] = si1[i];
    si13[2*i+1] = si3[i];
    si13e[2*i] = si1e[i];
    si13e[2*i+1] = si3e[i];
    si13o[2*i] = si1o[i];
    si13o[2*i+1] = si3o[i];
#endif // SHIFT13
  }
  //TRACE_ALL;
}
