typedef struct {
  int adj;
  int wantAnorm2;
  double Anorm2;
} LinopInfo;
#define LINOP_INFO_DEFAULT {0,0,0}

#define real float
#define PR(x) x ## F
#include "solversP.h"
#undef real
#undef PR

#define real double
#define PR(x) x ## D
#include "solversP.h"
#undef real
#undef PR
