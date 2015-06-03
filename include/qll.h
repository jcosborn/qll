#include "layout.h"
#include "thread.h"
#include "solvers.h"

void *myalloc(size_t size);

typedef struct {
  double one_link;
  double three_staple;
  double five_staple;
  double seven_staple;
  double lepage;
  double naik;
} Fat7Coeffs;
#define FAT7COEFFS_ZERO (Fat7Coeffs){0,0,0,0,0,0}

#define real float
#define PR(x) x ## F
#include "qllP.h"
#undef real
#undef PR

#define real double
#define PR(x) x ## D
#include "qllP.h"
#undef real
#undef PR

#define beginParallel BEGIN_PARALLEL
#define beginParallelSubset(s) beginParallel { SUBSET_DEFS2; setOffSub2(s);
#define endParallel END_PARALLEL
#define endParallelSubset } endParallel
#define singleThread _Pragma("omp single")

#define Field P(Field)
#define fieldNew P(fieldNew)
#define fieldFree P(fieldFree)
#define fieldSet P(fieldSet)
#define fieldZero P(fieldZero)
#define fieldSumElem P(fieldSumElem)
#define fieldNorm2 P(fieldNorm2)

#define Transporter P(Transporter)
#define transporterNew P(transporterNew)
#define transporterFree P(transporterFree)
#define transporterStartT P(transporterStartT)
#define transporterLocalT P(transporterLocalT)
#define transporterWaitT P(transporterWaitT)
#define transporterDoT P(transporterDoT)
#define transporterDo P(transporterDo)

#define gaugeDoublestore P(gaugeDoublestore)
#define gaugePlaq P(gaugePlaq)
#define smearFat7 P(smearFat7)

#define F_eq_r(ff,r) X_eq_r((ff)->f, (r), (ff)->nelem, ti0, ti1)
#define lnorm2_F(ff) norm2_X((ff)->f, (ff)->nelem, ti0, ti1)
#define sumDoubleArray(x,n) gsum(x,n,TC)
