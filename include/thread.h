#ifdef _OPENMP
#include <omp.h>
#define OMP_PARALLEL _Pragma("omp parallel")
#define OMP_NUMTHREADS omp_get_num_threads()
#define OMP_THREADNUM  omp_get_thread_num()
#define OMP_BARRIER _Pragma("omp barrier")
#else
#define OMP_PARALLEL
#define OMP_NUMTHREADS 1
#define OMP_THREADNUM  0
#define OMP_BARRIER
#endif
#define NUMTHREADS OMP_NUMTHREADS
#define THREADNUM OMP_THREADNUM
#define TBARRIER OMP_BARRIER
#ifndef TIBDEFAULT
#define TIBDEFAULT 0
#endif
#define MAXTHREADS 512
#define BEGIN_PARALLEL		    \
  {				    \
  double *_tsum[MAXTHREADS];	    \
  int volatile _tbar[MAXTHREADS];   \
  OMP_PARALLEL			    \
  {				    \
  DEFINE_THREAD_CTX;		    \
  _tbar[TC0.tid] = 0;		    \
  OMP_BARRIER;
#define END_PARALLEL }}

typedef struct {
  double **sum;
  int volatile *tbar;
  int tid;
  int nid;
  int tib;
  int n0, n1, d0;
} ThreadCtx;

#define DEFINE_THREAD_CTX			\
  ThreadCtx TC0, *TC;				\
  TC = &TC0;					\
  TC0.tid = OMP_THREADNUM;			\
  TC0.nid = OMP_NUMTHREADS;			\
  TC0.tbar = _tbar;				\
  TC0.sum = _tsum;				\
  SET_TIB(TIBDEFAULT)

#define SET_TIB(n)				\
  TC->tib = (n);				\
  TC->n0 = 10*TC->tid - TC->tib;		\
  TC->n1 = 10*(TC->tid+1) - TC->tib;		\
  if(TC->n0<0) TC->n0 = 0;			\
  TC->d0 = 10*TC->nid - TC->tib;		\

#define SUBSET_DEFS int ti0, ti1, off, len
#define SUBSET_DEFS2 int ti0, ti1
#define SET_SUBSET(s)
#define TI0(n) ((TC->tid*(n))/TC->nid)
#define TI1(n) (((TC->tid+1)*(n))/TC->nid)
//#define TI0C(n) TI0(n)
//#define TI1C(n) TI1(n)
//#define TI0C(n) (tid==0?0:(((tid-1)*(n))/(nid-1)))
//#define TI1C(n) (((tid)*(n))/(nid-1))
//#define TI0C(n) (tid<4?0:(((tid-4)*(n))/(nid-4)))
//#define TI1C(n) (((tid-3)*(n))/(nid-4))
//#define TI0C(n) (tid==0?0:(((2*tid-1)*(n))/(2*nid-1)))
//#define TI1C(n) (((2*tid+1)*(n))/(2*nid-1))
#define TI0C(n) ((TC->n0*(n))/TC->d0)
#define TI1C(n) ((TC->n1*(n))/TC->d0)
#define setOff(n) ti0=TI0(n); ti1=TI1(n); off=ti0; len=ti1-ti0
#define setOffC(n) ti0=TI0C(n); ti1=TI1C(n); off=ti0; len=ti1-ti0
//#define setOffS(n) ti0=TI0(n); ti1=TI1(n); off=ti0; len=ti1
#define setOffSub(s) setOff((s).endOuter-(s).beginOuter); ti0+=(s).beginOuter; ti1+=(s).beginOuter; off+=(s).beginOuter
#define setOffSubC(s) setOffC((s).endOuter-(s).beginOuter); ti0+=(s).beginOuter; ti1+=(s).beginOuter; off+=(s).beginOuter

#define setOff2(n) ti0=TI0(n); ti1=TI1(n)
#define setOffC2(n) ti0=TI0C(n); ti1=TI1C(n)
#define setOffSub2(s) setOff2((s).endOuter-(s).beginOuter); ti0+=(s).beginOuter; ti1+=(s).beginOuter
#define setOffSubC2(s) setOffC2((s).endOuter-(s).beginOuter); ti0+=(s).beginOuter; ti1+=(s).beginOuter


#ifdef __powerpc__
//#define SYNC __sync()
#define SYNC __eieio()
#else
#define SYNC
#endif

#define T0WAIT OMP_BARRIER

#if 0
#define T0WAIT					\
  if(tid==0) {					\
    tbar[0]++;					\
    int tbar0 = tbar[0];			\
    for(int b=1; b<nid; b++)			\
      while(tbar[b]<tbar0);			\
  } else {					\
    tbar[tid]++;				\
  }
#endif

#if 0
#define T0WAIT					\
  {						\
    int _tb1 = tbar[tid] + 1;			\
    int _b = 1;					\
    while((tid&_b)==0) {			\
      int _t1 = tid + _b;			\
      if(_t1>=nid) break;			\
      while(tbar[_t1]<_tb1);			\
      _b *= 2;					\
    }						\
    tbar[tid] = _tb1;				\
  }
#endif

#if 0
#define T0WAIT t0wait(tid, nid)
void static inline
t0wait(int tid, int nid)
{
  int _tb1 = tbar[tid] + 1;
  int _b = 1;
  while((tid&_b)==0) {
    int _t1 = tid + _b;
    if(_t1>=nid) break;
    while(tbar[_t1]<_tb1);
    _b *= 2;
  }
  tbar[tid] = _tb1;
}
#endif

//#define TWAIT0 OMP_BARRIER
#define TWAIT0 twait0(TC)

#if 0
#define TWAIT0					\
  if(tid==0) {					\
    tbar[0]++;					\
  } else {					\
    tbar[tid]++;				\
    int tbar0 = tbar[tid];			\
    while(tbar[0]<tbar0);			\
  }
#endif

void static inline
twait0(ThreadCtx *TC)
{
  SYNC;
  int volatile *tbar = TC->tbar;
  int tid = TC->tid;
  if(tid==0) {
    tbar[0]++;
  } else {
    tbar[tid]++;
    int tbar0 = tbar[tid];
    while(tbar[0]<tbar0);
  }
  SYNC;
}
