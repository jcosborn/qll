#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include <omp.h>
#include "common.h"

#define VpeqMtimesVXp(x,m,i,p,y) P(V_vpeq_M_times_xVp)(x,m,i,p,y,ti0,ti1)
#define FpeqMtimesFXp(x,m,i,p,y,n,c,t) P(X_peq_M_times_xXp)(x,m,i,p,y,n,c,t,ti0,ti1)
#define FXpeqMXtimesFXb(i0,x,m,b,i1,y1,i2,y2,n,c,t) P(xX_peq_xM_times_yXzXb)(i0,x,m,b,i1,y1,i2,y2,n,c,t,ti0,ti1)
#define VXpeqMXtimesVXb(i0,x,m,b,i1,y1,i2,y2) P(xV_vpeq_xM_times_yVzVb)(i0,x,m,b,i1,y1,i2,y2,ti0,ti1)

// U(n,n) F(o,n,m)
void
transporterNew(Transporter *t, real *u, ShiftIndices *si,
	       Subset *sub, int nelem, int nc)
{
  t->u = u;
  t->si = *si;
  t->nelem = nelem;
  t->nc = nc;
  prepareShiftBuf(&(t->sb), si, nelem*sizeof(real));
  t->sub = *sub;
  t->backward = 0;
}

void
transporterFree(Transporter *t)
{
  freeShiftBuf(&(t->sb));
}

// pTransStart(pt, dest, source, c) (dest+c*pt(src))
void
transporterStartT(Transporter *t, real *dest, real *src, int td, ThreadCtx *TC)
{
  SUBSET_DEFS;
  int tid = TC->tid;

  //DSTAMP(t0);
  if(tid==0) {
    if(t->si.nRecvRanks>0) {
      startRecvBuf(&(t->sb));
    }
  }

  //DSTAMP(t1);
  //double tp=0,tbs=0,ts=0;
  if(t->si.nSendRanks>0) {
    //DSTAMP(tt0);
    //setOff(t->si.nSendSites);
    setOffC(t->si.nSendSites);
    P(X_veq_pack_xX)(t->sb.sbuf, t->si.pack, t->si.sendSites, src,
		     t->nelem, off, len);
    //DSTAMP(tt1);
    T0WAIT;
    //DSTAMP(tt2);
    if(tid==0) {
      startSendBuf(&(t->sb));
      //printf("vvs: %i\tpack: %i\n", s->vvs, s->pack);
      //for(int i=0; i<8; i++) {
      //printf("idxSend[%i]: %i\n", i, s->idxSend[i]);
      //}
    }
    //DSTAMP(tt3);
    //tp += tt1-tt0;
    //tbs += tt2-tt1;
    //ts += tt3-tt2;
  }
  t->dest = dest;
  t->src = src;
  t->td = td;
  //t->scale = scale;
}

void
transporterLocalT(Transporter *t, ThreadCtx *TC)
{
  SUBSET_DEFS2;
  //DSTAMP(t2);
  setOffSubC2(t->sub);
  if(t->td>0 && t->u) {
    //printf("%i: %i %i\n", TC->tid, ti0, ti1);
    //printf("%i: %i %i %i\n", TC->tid, t->nelem, t->nc, t->td);
    //VpeqMtimesVXp(t->dest, t->u, t->si.pidx, t->si.perm, t->src);
    //printf("perm: %i\n", t->si.perm);
    FpeqMtimesFXp(t->dest, t->u, t->si.pidx, t->si.perm, t->src,
		  t->nelem, t->nc, t->td);
  } else {
    P(X_eq_xXp)(t->dest, t->si.pidx, t->si.perm, t->src, t->nelem, ti0, ti1);
  }
  //DSTAMP(t21);
}

void
transporterWaitT(Transporter *t, ThreadCtx *TC)
{
  SUBSET_DEFS2;
  int tid = TC->tid;
  //DSTAMP(t3);
  //double tws=0,twr=0,tbr=0,tb=0,tdr=0;
  if(t->si.nRecvDests>0) {
    //DSTAMP(tt0);
    //double tt1;
    if(t->si.nRecvRanks>0) {
      if(tid==0) {
	waitRecvBuf(&(t->sb));
      }
      //STAMP(tt1);
      TWAIT0;
      //} else {
      //STAMP(tt1);
    }
    //DSTAMP(tt2);
    if(t->sb.offr[tid]<0) {
      setOffSubC2(t->sub);
      int i0=0; while(i0<t->si.nRecvDests && t->si.recvDests[i0]<ti0) i0++;
      int i1=i0; while(i1<t->si.nRecvDests && t->si.recvDests[i1]<ti1) i1++;
      t->sb.offr[tid] = i0;
      t->sb.lenr[tid] = i1;
      //printf0("%i %i: i0: %i\ti1: %i\trd0: %i\trd1: %i\n", i, tid, i0, i1, t->si.recvDests[i0], (i1-1)<0?0:t->si.recvDests[i1-1]);
    }
    ti0 = t->sb.offr[tid];
    ti1 = t->sb.lenr[tid];
    if(t->td>0 && t->u) {
      //VXpeqMXtimesVXb(t->si.recvDests, t->dest, t->u, t->si.blend,
      //t->si.recvLocalSrcs, t->src, t->si.recvRemoteSrcs,
      //t->sb.rbuf);
      FXpeqMXtimesFXb(t->si.recvDests, t->dest, t->u, t->si.blend,
		      t->si.recvLocalSrcs, t->src, t->si.recvRemoteSrcs,
		      t->sb.rbuf, t->nelem, t->nc, t->td);
    } else {
      P(xX_eq_blend_xXxX)(t->si.recvDests, t->dest, t->si.blend,
			  t->si.recvLocalSrcs, t->src, t->si.recvRemoteSrcs,
			  t->sb.rbuf, t->nelem, ti0, ti1);
    }
    //DSTAMP(tt3);
    if(t->si.nRecvRanks>0) {
      if(tid==0) {
	doneRecvBuf(&(t->sb));
      }
    }
    //DSTAMP(tt4);
    //twr += tt1-tt0;
    //tbr += tt2-tt1;
    //tb += tt3-tt2;
    //tdr += tt4-tt3;
  }
  if(t->si.nSendRanks>0) {
    //DSTAMP(tt0);
    if(tid==0) {
      waitSendBuf(&(t->sb));
    }
    //DSTAMP(tt1);
    //tws += tt1-tt0;
  }

#if 0
  //DSTAMP(t4);
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
#endif
}

void
transporterDoT(Transporter *t, real *dest, real *src, int td, ThreadCtx *TC)
{
  transporterStartT(t, dest, src, td, TC);
  transporterLocalT(t, TC);
  transporterWaitT(t, TC);
}

void
transporterDo(Transporter *t, real *dest, real *src, int td)
{
  BEGIN_PARALLEL
    {
      transporterDoT(t, dest, src, td, TC);
    }
  END_PARALLEL;
}
