#include <stdlib.h>
#include <qmp.h>
#include "vla.h"
#include "qll.h"
//#include "thread.h"
#include "rdtsc.h"
//#include "solvers.h"

#define NAIK
#define SHIFT13

#define if0 if(myrank==0)
#define printf0(...) if(myrank==0) { printf(__VA_ARGS__); fflush(stdout); }
#define TRACE do{if0{printf("%s %s %i\n", __FILE__, __func__, __LINE__);fflush(stdout);}} while(0)
#define TRACET do{if0{printf("%i: %s %s %i\n", TC->tid, __FILE__, __func__, __LINE__);fflush(stdout);}} while(0)
#define TRACE_ALL do{printf("%i %s %s %i\n", myrank, __FILE__, __func__, __LINE__);fflush(stdout);} while(0)
#define BEGIN_FATAL (void)0
#define PRINT_FATAL(...) printf0(__VA_ARGS__)
#define END_FATAL exit(-1)

#define PRINTV(s,f,v,n) do { printf(s);			\
    for(int _i=0; _i<n; _i++) printf(" "f, (v)[_i]);	\
    printf("\n"); } while(0)
#define PRINTV0(s,f,v,n) do { if0{ PRINTV(s,f,v,n); }} while(0)

#define ARRAY_CREATE(tp,nm) int CAT(n,nm)=0, CAT(nm,len)=0; tp *nm=NULL
#define ARRAY_GROW(tp,nm,ln)			\
  if(CAT(n,nm)+(ln)>CAT(nm,len)) {		\
    CAT(nm,len) += (ln);			\
    nm = realloc(nm, CAT(nm,len)*sizeof(tp));	\
  }
#define ARRAY_APPEND(tp,nm,vl)			\
  if(CAT(n,nm)==CAT(nm,len)) {			\
    CAT(nm,len) *= 2;				\
    if(CAT(nm,len)==0) CAT(nm,len) = 16;	\
    nm = realloc(nm, CAT(nm,len)*sizeof(tp));	\
  }						\
  nm[CAT(n,nm)] = vl;				\
  CAT(n,nm)++
#define ARRAY_COPY(tp,ds,nm) \
  do { for(int _i=0; _i<CAT(n,nm); _i++) (ds)[_i] = (nm)[_i]; } while(0)
#define ARRAY_CLONE(tp,ds,nm) \
  do { (ds) = myalloc(CAT(n,nm)*sizeof(tp)); ARRAY_COPY(tp,ds,nm); } while(0)

extern int myrank;

typedef void GatherMap(int *srcRank, int *srcIdx, int dstRank, int *dstIdx, void *args);
void makeGather(GatherIndices *gi, GatherMap *map, void *args,
		int nSrcRanks, int nDstRanks, int myndi, int myrank);
void makeGatherFromGD(GatherIndices *gi, GatherDescription *gd);
void P(shift)(real *x, real *y, int nelem, int dir, int len);
ShiftIndices *P(getStaggeredSI)(int naik, char *sub);
ShiftIndices *P(getShift)(int dir, int len, char *sub);
void
P(createStaggeredLinksFb)(P(StaggeredLinks) *sl, Layout *l,
			  real *u[], real *u3[], char *sub);
void P(dslash)(real *x, real *u[8], real *u3[8], real mass, real *y, char *sub);
//void dslashT(real *x, real *u[8], real *u3[8], real mass, real *y, int nid, int tid);
void P(dslashFb)(real *x, real *u[8], real *u3[8], real mass, real *y, char *sub);
void P(dslashResetTimer)(void);
void P(dslashPrintTimer)(double scale);

void P(wilsonDslashResetTimer)(void);
void P(wilsonDslashPrintTimer)(double scale);

void qopStagDslashInit(int *argc, char ***argv);
void qopStagDslashFini(void);
void qopStagDslashSetup(Layout *l);
void qopShift(Layout *l, real *xx, real *yy, int nelem, int dir, int len);
void qopTransportM(Layout *l, real *xx, real *mm, real *yy, int dir, int len);
void qopSmear(Field *fl[], Field *ll[], Fat7Coeffs *coef, Field *u[], Field *ul[]);
void qopDslash(Layout *l, real *x, real *u[8], real *u3[8], real mass, real *y, char *sub);
void qopSolve(Layout *l, real *x, real *u[8], real *u3[8], real mass, real *y, double rsq, char *sub);
void
qopSolveMulti(Layout *l, real *x[], real *u[8], real *u3[8],
	      double masses[], real *y, int nmasses, double rsq, char *sub);

void qopWilsonDslashInit(int *argc, char ***argv);
void qopWilsonDslashFini(void);
void qopWilsonDslashSetup(Layout *l);
void qopSproj(Layout *l, real *x, real *y, int dir, int sign);
void qopSrecon(Layout *l, real *x, real *y, int dir, int sign);
void qopWilsonDslash(Layout *l, real *x, real *u[8], real mass, int sign, real *y, char *sub);
void qopWilsonSolve(Layout *l, real *x, real *u[8], real mass, real *y, double rsq, char *sub);
void qopWilsonSolveMulti(Layout *l, real *x[], real *u[8], double masses[],
			 real *y, int nmasses, double rsq, char *sub);

void makeShift(ShiftIndices *si, Layout *l, int *disp);
void makeShiftSub(ShiftIndices *si, Layout *l, int *disp, char *sub);
void makeShiftMulti(ShiftIndices *si[], Layout *l, int *disp[], int ndisp);
void makeShiftMultiSub(ShiftIndices *si[], Layout *l, int *disp[],
		       char *subs[], int ndisp);
void startSend(void *buf, int size, ShiftIndices *si);
void startRecv(void *buf, int size, ShiftIndices *si);
void waitSend(ShiftIndices *si);
void waitRecv(ShiftIndices *si);
void doneRecvBuf2(ShiftIndices *si);
void freeSend(ShiftIndices *si);
void freeRecv(ShiftIndices *si);

void prepareShiftBuf(ShiftBuf *sb, ShiftIndices *si, int esize);
void startSendBuf(ShiftBuf *sb);
void startRecvBuf(ShiftBuf *sb);
void waitSendBuf(ShiftBuf *sb);
void waitRecvBuf(ShiftBuf *sb);
void doneRecvBuf(ShiftBuf *sb);
void freeShiftBuf(ShiftBuf *sb);
void prepareShiftBufs(ShiftBuf *sb[], ShiftIndices *si[], int n, int esize);
void freeShiftBufs(ShiftBuf *sb[], int n);

void gsum(double *x, int n, ThreadCtx *TC);
